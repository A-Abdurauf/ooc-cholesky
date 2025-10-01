#include <sys/time.h>

#include <vector>
#include <mutex>

#include "common.h"
double get_current_time() {
  struct timeval time_val;

  struct timezone time_zone;

  gettimeofday(&time_val, &time_zone);

  return (double)(time_val.tv_sec) + (double)(time_val.tv_usec) / 1000000.0;
}

struct EventLog {
  int threadId;
  int deviceId;
  int task;
  double startTime;
  double endTime;
};

std::mutex mutex;

std::vector<EventLog> &get_event_log() {
  static std::vector<EventLog> event_log;
  return event_log;
}

void log_event(int rank, int device, int task, double start, double end) {
  EventLog eventLog{.threadId=rank, .deviceId=device, .task=task, .startTime=start, .endTime=end};
  std::lock_guard lock{mutex};
  get_event_log().push_back(eventLog);
}
