# FlashWeave (Tackmann, Matias Rodrigues & von Mering 2019, Cell Systems 9:286) on each layer: heterogeneous + sensitive mode
# (FlashWeaveHE-S: zero-adjusted partial-correlation tests on samples where both nodes are present, clr-adapt normalisation),
# max_k = 3, alpha = 0.01, with the meta-variables of netcommon.meta_vars as extra nodes, so that associations explained by a
# shared driver (winter, pool, medium, ice layer, day, under-ice light) are removed by conditioning.
# Usage: JULIA_DEPOT_PATH=... julia -p <workers> run_flashweave.jl <network dir>
# Reads cache/{layer}_X.tsv and cache/{layer}_meta.tsv; writes out/flashweave_{layer}.edgelist (all edges, meta included).
using Distributed
@everywhere using FlashWeave
net = length(ARGS) > 0 ? ARGS[1] : "."
layers = [replace(f, "_X.tsv" => "") for f in readdir(joinpath(net, "cache")) if endswith(f, "_X.tsv")]   # every layer export.py wrote
for layer in layers
    t = time()
    res = learn_network(joinpath(net, "cache", "$(layer)_X.tsv"), joinpath(net, "cache", "$(layer)_meta.tsv");
                        sensitive=true, heterogeneous=true, max_k=3, alpha=0.01, n_obs_min=-1, verbose=true)
    save_network(joinpath(net, "out", "flashweave_$(layer).edgelist"), res)
    println("$(layer): FlashWeave done in $(round(time() - t, digits=1)) s")
end
