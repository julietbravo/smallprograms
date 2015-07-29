# Julia implementation of diffusion

function diff(at, a, visc, dxidxi, dyidyi, dzidzi, itot, jtot, ktot)
    ii = 1
    jj = itot
    kk = itot*jtot

    @inbounds for k = 2:ktot-1
        for j = 2:jtot-1
            @simd for i = 2:itot-1
                ijk = i + (j-1)*jj + (k-1)*kk
                at[ijk] += visc * (
                        + ( (a[ijk+ii] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-ii]) ) * dxidxi 
                        + ( (a[ijk+jj] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
                        + ( (a[ijk+kk] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-kk]) ) * dzidzi
                        )
            end
        end
    end
end

function init(a, at, ncells)
    for i = 1:ncells
        a[i] = (i-1)^2 / i^2
        at[i] = 0
    end
end 

# Start benchmark
nloop  = 1000
itot   = 128
jtot   = 128
ktot   = 128
ncells = itot*jtot*ktot

a  = zeros(ncells)
at = zeros(ncells)

init(a, at, ncells)

# Check results
diff(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot)
println("at=",at[itot*jtot+itot+itot/2+1]);

# Time performance 
tic()

for i = 1:nloop
    diff(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot)
end

elapsed = toq()

println("time/iter =", elapsed/nloop, " s (", nloop, " iters)")
