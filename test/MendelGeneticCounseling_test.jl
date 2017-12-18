using MendelBase                        # From package MendelBase.
using SearchSetup                       # From package Search.
using DataFrames                        # From package DataFrames.
using MendelGeneticCounseling

@testset "penetrance_genetic_counseling" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "genetic counseling 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    # first test "genetic counseling 1" files
    multi_genotype = zeros(Int, 2, 1, 3) #need an Array{Int64,2} element
    multi_genotype[:, 1, 1] = [1; 2]
    multi_genotype[:, 1, 2] = [1; 1]
    multi_genotype[:, 1, 3] = [2; 1]

    par = parameter.par #this is [0.0]

    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)


    # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    # Here there are 2 pedigrees: first one is from 1 to 74, second one is fron 75 to 142
    n = 1
    start = instruction.extra[n][1]
    finish = instruction.extra[n][2]
    i = instruction.extra[n][3]

    for j = 1:3
        number = MendelGeneticCounseling.penetrance_genetic_counseling(person, locus, 
            multi_genotype[:, :, j], par, keyword, start, finish, i)
        @test number == 1.0
    end

    # for n = 1:73
    #     # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    #     # Here there are 2 pedigrees: first one is from 1 to 74, second one is fron 75 to 142
    #     # For some reason if n ≥ 74 the code will have bound array problems?
    #     start = instruction.extra[n][1]
    #     finish = instruction.extra[n][2]
    #     i = instruction.extra[n][3]
    #     if i == 0 continue end #Do not call penetrance_genetic_counseling if i is 0


    #     for j = 1:3
    #         number = MendelGeneticCounseling.penetrance_genetic_counseling(person, locus, 
    #             multi_genotype[:, :, j], par, keyword, start, finish, i)
    #         # @test number == 1.0
    #     end
    # end
end

@testset "prior_genetic_counseling" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "genetic counseling 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    # first test "genetic counseling 1" files
    multi_genotype = zeros(Int, 2, 1, 3) #need an Array{Int64,2} element
    multi_genotype[:, 1, 1] = [1; 2]
    multi_genotype[:, 1, 2] = [1; 1]
    multi_genotype[:, 1, 3] = [2; 1]

    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    par = parameter.par #this is [0.0]

    # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    # Here there are 2 pedigrees: first one is from 1 to 74, second one is fron 75 to 142
    n = 1
    start = instruction.extra[n][1]
    finish = instruction.extra[n][2]
    i = instruction.extra[n][3]

    prob = MendelGeneticCounseling.prior_genetic_counseling(person, locus, 
        multi_genotype[:, :, 3], par, keyword, start, finish, i)
end

@testset "transmission_genetic_counseling" begin
    1 == 1
end

@testset "basic and wrappers" begin
    1 == 1
end
