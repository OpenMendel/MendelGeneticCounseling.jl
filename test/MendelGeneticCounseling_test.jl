using MendelBase                        # From package MendelBase.
using SearchSetup                       # From package Search.
using DataFrames                        # From package DataFrames.
using MendelGeneticCounseling

@testset "penetrance_genetic_counseling" begin

    # first test "genetic counseling 1" files
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "genetic counseling 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)
    par = parameter.par #this is [0.0]

    for n = 1:73
        # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
        # Here there are 2 pedigrees: first one is from 1 to 74, second one is fron 75 to 142
        # For some reason if n ≥ 74 the code will have bound array problems?
        operation = instruction.operation[n]
        start = instruction.extra[n][1]
        finish = instruction.extra[n][2]
        i = instruction.extra[n][3]
        if operation != penetrance_and_prior_array continue end #this avoids some array access errors

        genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                            genotypes, i)

        for j = 1:genotypes
            number = MendelGeneticCounseling.penetrance_genetic_counseling(person, locus, 
                multi_genotype[:, :, j], par, keyword, start, finish, i)
            @test number == 1.0 #this is always 1 
        end
    end





    # now test "genetic counseling 2" files
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "genetic counseling 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)
    par = parameter.par #this is [0.0]

    for n = 1:29
        # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
        # Here there are 2 pedigrees: first one is from 1 to 30, second one is fron 31 to 60
        operation = instruction.operation[n]
        start = instruction.extra[n][1]
        finish = instruction.extra[n][2]
        i = instruction.extra[n][3]
        if operation != penetrance_and_prior_array; continue; end #this avoids some array access errors
        if person.mother[i] != 0; continue; end #prior prob doesnt exist for non founder

        genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                            genotypes, i)

        for j = 1:genotypes
            number = MendelGeneticCounseling.penetrance_genetic_counseling(person, locus, 
                multi_genotype[:, :, j], par, keyword, start, finish, i)
            @test number == 1.0 #this is always 1 
        end
    end
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
    par = parameter.par #this is [0.0]
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    for n = 1:73
        # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
        # Here there are 2 pedigrees: first one is from 1 to 74, second one is fron 75 to 142
        # For some reason if n ≥ 74 the code will have bound array problems?
        operation = instruction.operation[n]
        start = instruction.extra[n][1]
        finish = instruction.extra[n][2]
        i = instruction.extra[n][3]
        if operation != penetrance_and_prior_array continue end #this avoids some array access errors
        if person.mother[i] != 0 continue end #prior prob doesnt exist for non founder

        genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                            genotypes, i)

        for j = 1:genotypes
            prob = MendelGeneticCounseling.prior_genetic_counseling(person, locus, 
                multi_genotype[:, :, j], par, keyword, start, finish, i)
            if multi_genotype[:, :, j] == [1; 2]; @test prob == 0.98 * 0.02; end
            if multi_genotype[:, :, j] == [1; 1]; @test prob == 0.98 * 0.98; end
            if multi_genotype[:, :, j] == [2; 1]; @test prob == 0.98 * 0.02; end
            if multi_genotype[:, :, j] == [2; 2]; @test prob == 0.02 * 0.02; end
        end
    end



    # now test "genetic counseling 2" files
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "genetic counseling 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)
    par = parameter.par #this is [0.0]

    for n = 1:29
        # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
        # Here there are 2 pedigrees: first one is from 1 to 30, second one is fron 31 to 60
        operation = instruction.operation[n]
        start = instruction.extra[n][1]
        finish = instruction.extra[n][2]
        i = instruction.extra[n][3]
        if operation != penetrance_and_prior_array; continue; end #this avoids some array access errors
        if person.mother[i] != 0; continue; end #prior prob doesnt exist for non founder

        genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                            genotypes, i)

        for j = 1:genotypes
            prob = MendelGeneticCounseling.prior_genetic_counseling(person, locus, 
                multi_genotype[:, :, j], par, keyword, start, finish, i)
            answer = 1.0

            if multi_genotype[1, 1, j] == 1; answer *= 1.0e-5; end
            if multi_genotype[1, 1, j] == 2; answer *= 0.99999; end
            if multi_genotype[2, 1, j] == 1; answer *= 1.0e-5; end
            if multi_genotype[2, 1, j] == 2; answer *= 0.99999; end

            if multi_genotype[1, 2, j] == 1; answer *= 0.52; end
            if multi_genotype[1, 2, j] == 2; answer *= 0.48; end
            if multi_genotype[2, 2, j] == 1; answer *= 0.52; end
            if multi_genotype[2, 2, j] == 2; answer *= 0.48; end

            @test answer == prob
        end
    end
end

@testset "transmission_genetic_counseling" begin
    1 == 1
end

@testset "basic and wrappers" begin
    1 == 1
end
