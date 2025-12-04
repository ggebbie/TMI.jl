# ⚠️ This version leaks in Julia 1.11+ but not in 1.10

function inner(x; z = nothing)
    return sum(abs2, x) + (z === nothing ? 0 : z)
end

# closure that binds a keyword argument
f = x -> inner(x; z = nothing)

x = rand(1000)

GC.gc()
mem_before = Sys.maxrss()

for i = 1:20_000
    f(x)
    if i % 5000 == 0
        GC.gc()
        @info "iter=$i  rss=$(Sys.maxrss() - mem_before) bytes leaked"
    end
end



# ✔ positional argument version — no leak on 1.10, 1.11, 1.12+

function inner2(x, z)
    return sum(abs2, x) + (z === nothing ? 0 : z)
end

struct Wrapper
    z
end

# callable struct — no keyword call, no closure capture
(w::Wrapper)(x) = inner2(x, w.z)

f2 = Wrapper(nothing)

x = rand(1000)

GC.gc()
mem_before = Sys.maxrss()

for i = 1:20_000
    f2(x)
    if i % 5000 == 0
        GC.gc()
        @info "iter=$i  rss=$(Sys.maxrss() - mem_before) bytes leaked"
    end
end
