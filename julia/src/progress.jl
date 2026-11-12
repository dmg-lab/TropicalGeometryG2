# Progress of long census steps: redrawn in a terminal, written as lines so that a log file stays readable.

"""
    TaskProgress(description, total; interval=300, live=stdout isa Base.TTY)

Report a loop of `total` steps on standard output: a bar, the count, the
elapsed time and an estimate of the time left, and a last line with the total
time. In a terminal (`live`) the bar is drawn at once and redrawn in place a
few times a second. Otherwise a new line is written every `interval` seconds:
ProgressMeter redraws its bar with carriage returns and escape codes, which a
Slurm log keeps as one ever-growing line. `advance!` is thread safe.
"""
struct TaskProgress
    description::String
    total::Int
    interval::Float64
    live::Bool
    started::Float64
    done::Threads.Atomic{Int}
    deadline::Threads.Atomic{Float64}
    finished::Threads.Atomic{Bool}
    output::ReentrantLock
    # Nanoseconds this process had spent collecting garbage when the loop started.
    collected::Int
end

"Seconds between two redraws of a live bar."
const LIVE_PROGRESS_INTERVAL = 0.25

function TaskProgress(description::AbstractString, total::Integer; interval::Real=300,
                      live::Bool=stdout isa Base.TTY)
    started = time()
    period = live ? LIVE_PROGRESS_INTERVAL : Float64(interval)
    progress = TaskProgress(String(description), Int(total), period, live, started,
                            Threads.Atomic{Int}(0), Threads.Atomic{Float64}(started + period),
                            Threads.Atomic{Bool}(false), ReentrantLock(), Base.gc_num().total_time)
    # The first step may be the longest one.
    live && total > 0 && _report(progress, 0, started)
    return progress
end

"Advance `progress` by `steps`, reporting when its interval has passed and finishing at its total."
function advance!(progress::TaskProgress, steps::Integer=1)
    done = Threads.atomic_add!(progress.done, Int(steps)) + Int(steps)
    done >= progress.total && return finish!(progress)
    now = time()
    deadline = progress.deadline[]
    # Of several threads past the deadline, only the one moving it reports.
    if now >= deadline && Threads.atomic_cas!(progress.deadline, deadline, now + progress.interval) === deadline
        _report(progress, done, now)
    end
    return nothing
end

"Print the last line of `progress` once, also when its loop stops early."
function finish!(progress::TaskProgress)
    lock(progress.output) do
        Threads.atomic_xchg!(progress.finished, true) && return
        progress.total == 0 && return
        done = progress.done[]
        elapsed = _duration(time() - progress.started)
        _draw(progress, string(progress.description, " ", _progress_bar(done, progress.total), "  ",
                               done >= progress.total ? "$(_count(progress.total)) in $elapsed" :
                               "stopped at $(_count(done))/$(_count(progress.total)) after $elapsed");
              last=true)
    end
    return nothing
end

"Print `message` as its own line during the loop of `progress`, keeping a live bar below it."
function note!(progress::TaskProgress, message)
    lock(progress.output) do
        live = progress.live && !progress.finished[]
        live && print("\r\e[K")
        println("  ", message)
        flush(stdout)
        live && _report(progress, progress.done[], time())
    end
    return nothing
end

function _report(progress::TaskProgress, done, now)
    elapsed = now - progress.started
    line = string(progress.description, " ", _progress_bar(done, progress.total), " ",
                  lpad(floor(Int, 100 * done / progress.total), 3), "%  ",
                  _count(done), "/", _count(progress.total), "  ", _duration(elapsed), " elapsed",
                  done > 0 ? ", ~$(_duration(elapsed * (progress.total - done) / done)) left" : "",
                  progress.live ? "" : _resident_memory(; live=true) * _group_memory() *
                                       _collection_share(progress, elapsed))
    # A report racing the last line must not be drawn after it.
    lock(() -> progress.finished[] || _draw(progress, line), progress.output)
    return nothing
end

"The share of the loop this process spent collecting garbage, when it is a tenth or more: its workers wait meanwhile."
function _collection_share(progress::TaskProgress, elapsed)
    share = (Base.gc_num().total_time - progress.collected) / 1e9 / max(elapsed, 1e-3)
    return share < 0.1 ? "" : ", $(round(Int, 100 * share))% of the time collecting garbage"
end

"Write `line`, over the previous one when `progress` is live; `last` ends a live bar."
function _draw(progress::TaskProgress, line; last::Bool=false)
    if progress.live
        print("\r  ", line, "\e[K", last ? "\n" : "")
    else
        println("  ", line)
    end
    flush(stdout)
    return nothing
end

"A function advancing a new `TaskProgress` by its argument, for code that does not know the type."
function progress_steps(description, total)
    progress = TaskProgress(description, total)
    return steps -> advance!(progress, steps)
end

"Run `f`, printing `description` before and its duration after: for long steps with no loop to count."
function timed(f, description::AbstractString)
    println("  ", description, "...")
    flush(stdout)
    started = time()
    result = f()
    println("  ", description, ": ", _duration(time() - started))
    flush(stdout)
    return result
end

"""
Run `command`, which must not throw when it fails (use `ignorestatus`), and
print every `interval` seconds how long it has run and its resident memory.
"""
function run_reporting(command, description::AbstractString; interval::Real=600)
    process = run(command; wait=false)
    started = time()
    while timedwait(() -> process_exited(process), interval; pollint=1.0) === :timed_out
        println("  ", description, " running for ", _duration(time() - started),
                _resident_memory(() -> getpid(process)), _group_memory())
        flush(stdout)
    end
    wait(process)
    return process
end

"""
Resident memory of a process; for this one, with `live`, also what its Julia
objects take: the rest is memory that was freed but not returned to the system,
or that of libraries.
"""
function _resident_memory(pid=getpid; live::Bool=false)
    try
        for line in eachline("/proc/$(pid())/status")
            startswith(line, "VmRSS:") &&
                return ", $(round(parse(Int, split(line)[2]) / 2^20; digits=1)) GB resident" *
                       (live ? " ($(_gigabytes(Base.gc_live_bytes())) live)" : "")
        end
    catch
        # The process has just exited, or this is not Linux.
    end
    return ""
end

"""
Memory used and allowed in the control group that limits this process, as a
Slurm job does: it counts the workers and gfan, and is what an out-of-memory
kill is decided on. The reclaimable part is the file cache not in use, which
counts as used until the kernel takes it back, before any kill.
"""
function _group_memory()
    try
        group = joinpath("/sys/fs/cgroup", lstrip(split(readline("/proc/self/cgroup"), "::")[2], '/'))
        while group != "/sys/fs/cgroup"
            limit = tryparse(Int, strip(read(joinpath(group, "memory.max"), String)))
            if !isnothing(limit)
                used = parse(Int, strip(read(joinpath(group, "memory.current"), String)))
                reclaimable = sum(parse(Int, split(line)[2]) for line in eachline(joinpath(group, "memory.stat"))
                                  if startswith(line, "inactive_file "); init=0)
                return ", $(_gigabytes(used)) of $(_gigabytes(limit)) GB in the job" *
                       " ($(_gigabytes(reclaimable)) reclaimable)"
            end
            group = dirname(group)
        end
    catch
        # No unified control group hierarchy, or not Linux.
    end
    return ""
end

_gigabytes(bytes) = round(bytes / 2^30; digits=1)

function _progress_bar(done, total; width=20)
    filled = total <= 0 ? width : clamp(floor(Int, width * done / total), 0, width)
    return "[" * "█"^filled * "░"^(width - filled) * "]"
end

_count(n) = replace(string(n), r"(?<=\d)(?=(\d{3})+$)" => ",")

function _duration(seconds)
    seconds < 10 && return "$(round(seconds; digits=1))s"
    s = round(Int, seconds)
    s < 60 && return "$(s)s"
    s < 3600 && return "$(s ÷ 60)m$(lpad(s % 60, 2, '0'))s"
    s < 86400 && return "$(s ÷ 3600)h$(lpad(s % 3600 ÷ 60, 2, '0'))m"
    return "$(s ÷ 86400)d$(lpad(s % 86400 ÷ 3600, 2, '0'))h"
end
