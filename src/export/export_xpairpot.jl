
function export_matrixpairpot(fname, IP)
    V1 = IP.components[1]
    V2 = IP.components[2]

    species = collect(string.(keys(V1.E0)))
    species_dict = Dict(zip(collect(0:length(species)-1), species))
    reversed_species_dict = Dict(zip(species, collect(0:length(species)-1)))

    
    data = Dict()    
    data["deltaSplineBins"] = 0.001

    elements = Vector(undef, length(species))
    E0 = zeros(3)
    
    for (index, element) in species_dict
        E0[index+1] = V1(Symbol(element))
        elements[index+1] = element
    end

    data["elements"] = elements
    data["E0"] = E0


    # need to make a dummy N=1 ACE with the same number of basis functions
    # type of V2 must be an Xpolypairpot

    pairdeg = length(V2.basis[1,1].J.A)
    rcut = 5.0
    println(pairdeg)
    println(rcut)

    dummy_ace = ACE1.Utils.rpi_basis(species = collect(keys(V1.E0)),
                            N = 1,
                            r0 = 1.0,
                            maxdeg = pairdeg,
                            rcut = rcut);

    
    embeddings, bonds = export_radial_basis(V3, species_dict)
    data["embeddings"] = embeddings
    data["bonds"] = bonds

    functions, lmax = export_ACE_functions(V3, species, reversed_species_dict)
    data["functions"] = functions
    data["lmax"] = lmax
end