#define _GNU_SOURCE

#include <assert.h>
#include <ctype.h>
#include <pthread.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW 50000.0
#define ALIENS_HIGH 150000.0

typedef struct {
  signal *sig;
  int filter_order;
  int num_bands;
  int num_processors;
  double bandwidth;
  double *band_power;
  int *next_band;
  pthread_mutex_t *next_band_lock;
  int thread_id;
} worker_args;

void usage()
{
  printf("usage: p_band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors\n");
}

double avg_power(double *data, int num)
{
  int i;
  double ss;

  ss = 0;
  for (i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double *data, int num)
{
  double m = data[0];
  int i;

  for (i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double *data, int num)
{
  double s = 0;
  int i;

  for (i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double *data, int num)
{
  int i;
  double dc = avg_of(data, num);

  printf("Removing DC component of %lf\n", dc);

  for (i = 0; i < num; i++) {
    data[i] -= dc;
  }
}

static double convolve_and_compute_power(int length, double input_signal[],
                                         int order, double coeffs[])
{
  int i, j;
  double sum_power = 0.0;

  for (i = 0; i < length; i++) {
    double output = 0.0;
    int max_j = i < order ? i : order;

    for (j = max_j; j >= 0; j--) {
      output += input_signal[i - j] * coeffs[j];
    }

    sum_power += output * output;
  }

  return sum_power / length;
}

static void bind_to_processor(int thread_id, int num_processors)
{
  cpu_set_t set;

  if (num_processors <= 0) {
    return;
  }

  CPU_ZERO(&set);
  CPU_SET(thread_id % num_processors, &set);
  sched_setaffinity(0, sizeof(set), &set);
}

static void *worker(void *arg)
{
  worker_args *args = (worker_args *)arg;
  double *filter_coeffs;

  bind_to_processor(args->thread_id, args->num_processors);

  filter_coeffs = (double *)malloc(sizeof(double) * (args->filter_order + 1));
  if (!filter_coeffs) {
    fprintf(stderr, "Out of memory\n");
    pthread_exit(NULL);
  }

  while (1) {
    int band;

    pthread_mutex_lock(args->next_band_lock);
    band = *args->next_band;
    (*args->next_band)++;
    pthread_mutex_unlock(args->next_band_lock);

    if (band >= args->num_bands) {
      break;
    }

    generate_band_pass(args->sig->Fs,
                       band * args->bandwidth + 0.0001,
                       (band + 1) * args->bandwidth - 0.0001,
                       args->filter_order,
                       filter_coeffs);
    hamming_window(args->filter_order, filter_coeffs);

    args->band_power[band] =
        convolve_and_compute_power(args->sig->num_samples,
                                   args->sig->data,
                                   args->filter_order,
                                   filter_coeffs);
  }

  free(filter_coeffs);
  pthread_exit(NULL);
}

int analyze_signal(signal *sig, int filter_order, int num_bands,
                   int num_threads, int num_processors, double *lb, double *ub)
{
  double Fc, bandwidth;
  double signal_power;
  double *band_power;
  pthread_t *threads;
  worker_args *args;
  pthread_mutex_t next_band_lock;
  int next_band = 0;
  int threads_to_start;
  double start, end;
  unsigned long long tstart, tend;
  resources rstart, rend, rdiff;
  int band;
  int rc;

  Fc = (sig->Fs) / 2;
  bandwidth = Fc / num_bands;

  remove_dc(sig->data, sig->num_samples);

  signal_power = avg_power(sig->data, sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  band_power = (double *)calloc(num_bands, sizeof(double));
  threads_to_start = num_threads < num_bands ? num_threads : num_bands;
  threads = (pthread_t *)malloc(sizeof(pthread_t) * threads_to_start);
  args = (worker_args *)malloc(sizeof(worker_args) * threads_to_start);

  if (!band_power || !threads || !args) {
    printf("Out of memory\n");
    free(band_power);
    free(threads);
    free(args);
    return 0;
  }

  pthread_mutex_init(&next_band_lock, NULL);

  get_resources(&rstart, THIS_PROCESS);
  start = get_seconds();
  tstart = get_cycle_count();

  for (band = 0; band < threads_to_start; band++) {
    args[band].sig = sig;
    args[band].filter_order = filter_order;
    args[band].num_bands = num_bands;
    args[band].num_processors = num_processors;
    args[band].bandwidth = bandwidth;
    args[band].band_power = band_power;
    args[band].next_band = &next_band;
    args[band].next_band_lock = &next_band_lock;
    args[band].thread_id = band;

    rc = pthread_create(&threads[band], NULL, worker, &args[band]);
    if (rc != 0) {
      fprintf(stderr, "Failed to start thread %d\n", band);
      exit(-1);
    }
  }

  for (band = 0; band < threads_to_start; band++) {
    rc = pthread_join(threads[band], NULL);
    if (rc != 0) {
      fprintf(stderr, "Failed to join thread %d\n", band);
      exit(-1);
    }
  }

  tend = get_cycle_count();
  end = get_seconds();
  get_resources(&rend, THIS_PROCESS);

  get_resources_diff(&rstart, &rend, &rdiff);

  double max_band_power = max_of(band_power, num_bands);
  double avg_band_power = avg_of(band_power, num_bands);
  int i;
  int wow = 0;

  *lb = *ub = -1;

  for (band = 0; band < num_bands; band++) {
    double band_low = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n"
         "User time        %lf seconds\n"
         "System time      %lf seconds\n"
         "Page faults      %ld\n"
         "Page swaps       %ld\n"
         "Blocks of I/O    %ld\n"
         "Signals caught   %ld\n"
         "Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\nNote that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  pthread_mutex_destroy(&next_band_lock);
  free(band_power);
  free(threads);
  free(args);

  return wow;
}

int main(int argc, char *argv[])
{
  signal *sig;
  double Fs;
  char sig_type;
  char *sig_file;
  int filter_order;
  int num_bands;
  int num_threads;
  int num_processors;
  double start, end;

  if (argc != 8) {
    usage();
    return -1;
  }

  sig_type = toupper(argv[1][0]);
  sig_file = argv[2];
  Fs = atof(argv[3]);
  filter_order = atoi(argv[4]);
  num_bands = atoi(argv[5]);
  num_threads = atoi(argv[6]);
  num_processors = atoi(argv[7]);

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);
  assert(num_threads > 0);
  assert(num_processors > 0);

  printf("type:     %s\n"
         "file:     %s\n"
         "Fs:       %lf Hz\n"
         "order:    %d\n"
         "bands:    %d\n"
         "threads:  %d\n"
         "procs:    %d\n",
         sig_type == 'T' ? "Text" : sig_type == 'B' ? "Binary" : sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE",
         sig_file,
         Fs,
         filter_order,
         num_bands,
         num_threads,
         num_processors);

  printf("Load or map file\n");

  switch (sig_type) {
  case 'T':
    sig = load_text_format_signal(sig_file);
    break;

  case 'B':
    sig = load_binary_format_signal(sig_file);
    break;

  case 'M':
    sig = map_binary_format_signal(sig_file);
    break;

  default:
    printf("Unknown signal type\n");
    return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  if (analyze_signal(sig, filter_order, num_bands, num_threads,
                     num_processors, &start, &end)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}
