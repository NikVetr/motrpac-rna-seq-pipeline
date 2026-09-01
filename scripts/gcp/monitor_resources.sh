#!/usr/bin/env bash
set -euo pipefail

interval_seconds="${MOTRPAC_MONITOR_INTERVAL_SECONDS:-15}"
max_samples="${MOTRPAC_MONITOR_MAX_SAMPLES:-0}"

case "${interval_seconds}:${max_samples}" in
    *[!0-9:]*|:*|*:)
        echo "Monitoring interval and sample count must be nonnegative integers" >&2
        exit 2
        ;;
esac

read_metric() {
    if [[ -r "$1" ]]; then
        awk 'NR == 1 { print $1 }' "$1"
    else
        printf 'NA\n'
    fi
}

cpu_usage_usec() {
    if [[ -r /sys/fs/cgroup/cpu.stat ]]; then
        awk '$1 == "usage_usec" { print $2 }' /sys/fs/cgroup/cpu.stat
    elif [[ -r /sys/fs/cgroup/cpuacct/cpuacct.usage ]]; then
        awk '{ printf "%.0f\n", $1 / 1000 }' \
            /sys/fs/cgroup/cpuacct/cpuacct.usage
    else
        printf 'NA\n'
    fi
}

memory_metric() {
    local v2_path="$1"
    local v1_path="$2"
    if [[ -r "$v2_path" ]]; then
        read_metric "$v2_path"
    else
        read_metric "$v1_path"
    fi
}

printf 'timestamp_utc\tepoch_s\tcpu_usage_usec\tmemory_current_bytes\tmemory_peak_bytes\tmemory_limit_bytes\thost_mem_available_kb\tdisk_used_kb\tdisk_available_kb\n'

sample_number=0
while true; do
    timestamp_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    epoch_s="$(date -u +%s)"
    cpu_usec="$(cpu_usage_usec)"
    memory_current="$(memory_metric \
        /sys/fs/cgroup/memory.current \
        /sys/fs/cgroup/memory/memory.usage_in_bytes)"
    memory_peak="$(memory_metric \
        /sys/fs/cgroup/memory.peak \
        /sys/fs/cgroup/memory/memory.max_usage_in_bytes)"
    memory_limit="$(memory_metric \
        /sys/fs/cgroup/memory.max \
        /sys/fs/cgroup/memory/memory.limit_in_bytes)"
    host_mem_available="$(awk \
        '$1 == "MemAvailable:" { print $2 }' /proc/meminfo)"
    read -r disk_used disk_available < <(
        df -Pk . | awk 'NR == 2 { print $3, $4 }'
    )

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$timestamp_utc" \
        "$epoch_s" \
        "$cpu_usec" \
        "$memory_current" \
        "$memory_peak" \
        "$memory_limit" \
        "$host_mem_available" \
        "$disk_used" \
        "$disk_available"

    sample_number=$((sample_number + 1))
    if (( max_samples > 0 && sample_number >= max_samples )); then
        break
    fi
    sleep "$interval_seconds"
done
