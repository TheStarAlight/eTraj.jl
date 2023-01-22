using .MolecularCalculators
using HDF5
using WignerD
using Base.Threads

"Sample provider which yields electron samples through WFAT formula, matching `IonRateMethod=:WFAT`"
struct WFATSampleProvider <: ElectronSampleProvider
    laser           ::Laser;
    target          ::Molecule;     # WFAT only supports [Molecule]
    tSamples        ::AbstractVector;
    ss_pdSamples    ::AbstractVector;
    ss_pzSamples    ::AbstractVector;
    wfat_intdata    ::Array;
    wfat_Ip         ::Real;
    wfat_μ          ::Vector;
    wfat_nξMax      ::Int;
    wfat_mMax       ::Int;
    wfat_lMax       ::Int;

    function WFATSampleProvider(;
                                laser               ::Laser,
                                target              ::Molecule,
                                sample_tSpan        ::Tuple{<:Real,<:Real},
                                sample_tSampleNum   ::Int,
                                ss_pdMax            ::Real,
                                ss_pdNum            ::Int,
                                ss_pzMax            ::Real,
                                ss_pzNum            ::Int,
                                wfat_intdata_path   ::String,
                                kwargs...   # kwargs are surplus params.
                                )
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "[WFATSampleProvider] Invalid time sample number $sample_tSampleNum."
        @assert (ss_pdNum>0 && ss_pzNum>0) "[WFATSampleProvider] Invalid pd/pz sample number $ss_pdNum/$ss_pzNum."
        # load WFAT IntData.
        wfat_intdata_file = nothing; wfat_intdata = nothing; molName = "";
        wfat_calcParams = nothing   # selectedOrbit_energy, μ_vec, nξMax, mMax, lMax
        if ! isfile(wfat_intdata_path)
            error("[WFATSampleProvider] WFAT Intdata \"$wfat_intdata_path\" doesn't exist.")
        end
        try
            wfat_intdata_file = h5open(wfat_intdata_path, "r")
            molName = read(wfat_intdata_file, "molName")
            wfat_intdata = read(wfat_intdata_file, "IntData")
            wfat_calcParams = read(wfat_intdata_file, "selectedOrbit_energy"), read(wfat_intdata_file, "μ_vec"), read(wfat_intdata_file, "nξMax"), read(wfat_intdata_file, "mMax"), read(wfat_intdata_file, "lMax")
        catch
            @error "[WFATSampleProvider] Encountered error when trying to read the WFAT Intdata."
            rethrow()
        finally
            if ! isnothing(wfat_intdata_file)
                close(wfat_intdata_file)
            end
        end
        if molName != target.name
            @warn "The molecule's name [$(target.name)] is different from that in the Intdata [$molName]. Check out if your data is correct."
        end
        # finish initialization.
        return new( laser, target,
                    range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                    range(-abs(ss_pdMax),abs(ss_pdMax);length=ss_pdNum), range(-abs(ss_pzMax),abs(ss_pzMax);length=ss_pzNum),
                    wfat_intdata, -1*wfat_calcParams[1], wfat_calcParams[2], wfat_calcParams[3], wfat_calcParams[4], wfat_calcParams[5]
                    )
    end
end

"Gets the total number of batches."
function batchNum(sp::WFATSampleProvider)
    return length(sp.tSamples)
end

"Gets the structure factor of the molecule target of a given channel under a specific Euler angle (α,β,γ)."
function getStructFactor(sp::WFATSampleProvider, nξ::Int, m::Int, α, β, γ)
    @assert nξ≥0 "[WFATSampleProvider] The nξ must be positive."
    sum = zero(ComplexF64)
    for l in abs(m):sp.wfat_lMax
        for m_ in -l:l
            sum += sp.wfat_intdata[nξ+1,m+sp.wfat_mMax+1,l+1,m_+l+1] * WignerD.wignerdjmn(l,m,m_,β) * exp(-1im*m_*γ)
        end
    end
    return real(sum)
end
