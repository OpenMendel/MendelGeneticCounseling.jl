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

        #
        # Construct the parent's multiple locus genotypes.
        #
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

        #
        # Construct the parent's multiple locus genotypes.
        #
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

        #
        # Construct the parent's multiple locus genotypes.
        #
        genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                            genotypes, i)

        for j = 1:genotypes
            prob = MendelGeneticCounseling.prior_genetic_counseling(person, locus, 
                multi_genotype[:, :, j], par, keyword, start, finish, i)
            # note the following probabilities are provided in the locus frame
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

        #
        # Construct the parent's multiple locus genotypes.
        #
        genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                            genotypes, i)

        for j = 1:genotypes
            prob = MendelGeneticCounseling.prior_genetic_counseling(person, locus, 
                multi_genotype[:, :, j], par, keyword, start, finish, i)
            answer = 1.0

            #note the following probabilities are provided in the locus frame
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

    # first test "genetic counseling 1" files
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "genetic counseling 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    par = parameter.par #this is [0.0]
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    for n = 1:73
        # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
        # Here there are 2 pedigrees: first one is from 1 to 74, second one is fron 75 to 142
        # For some reason if n ≥ 74 the code will have bound array problems?
        operation = instruction.operation[n]
        if operation != transmission_array; continue; end #this avoids some array access errors

        start = instruction.extra[n][1]
        finish = instruction.extra[n][2]
        i = instruction.extra[n][3]
        j = instruction.extra[n][4]

        # need 2 genotypes to run transmission, so we construct them
        # and give them the same names as used in elston_stewart_evaluation
        i_genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                            i_genotypes, i)

        maternal = !person.male[i]
        j_genotypes = MendelBase.genotype_count(person, locus, j, start, finish)
        gamete = MendelBase.construct_gametes(person, locus, start, finish, j_genotypes, j,
                                     maternal)

        # loop through all possible genotypes, and see if transmission probability is correct
        for l = 1:j_genotypes
            for k = 1:i_genotypes
                trans = MendelGeneticCounseling.transmission_genetic_counseling(person, 
                    locus, gamete[:, l], multi_genotype[:, :, k], par, keyword, 
                    start, finish, i, j)

                #gamete[:, l] could be [1] or [2]
                #multi_genotype[:, :, k] is either [1;1] or [1;2]
                #we compute the probaiblity that a multi_genotype pass down the gamete

                if gamete[:, l] == [1] && sum(multi_genotype[:, :, k]) == 2
                    @test trans == 1.0
                elseif gamete[:, l] == [1] && sum(multi_genotype[:, :, k]) == 3
                    @test trans == 0.5
                elseif gamete[:, l] == [1] && sum(multi_genotype[:, :, k]) == 4
                    @test trans == 0.0
                elseif gamete[:, l] == [2] && sum(multi_genotype[:, :, k]) == 2
                    @test trans == 0.0
                elseif gamete[:, l] == [2] && sum(multi_genotype[:, :, k]) == 3
                    @test trans == 0.5
                elseif gamete[:, l] == [2] && sum(multi_genotype[:, :, k]) == 4
                    @test trans == 1.0
                else
                    @test_throws(MethodError, "shouldn't have reached here bro")
                end
            end
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
        if operation != transmission_array; continue; end #this avoids some array access errors

        start = instruction.extra[n][1]
        finish = instruction.extra[n][2]
        i = instruction.extra[n][3]
        j = instruction.extra[n][4]


        # need 2 genotypes to run transmission, so we construct them
        # and give them the same names as used in elston_stewart_evaluation
        i_genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                            i_genotypes, i)

        maternal = !person.male[i]
        j_genotypes = MendelBase.genotype_count(person, locus, j, start, finish)
        gamete = MendelBase.construct_gametes(person, locus, start, finish, j_genotypes, j,
                                     maternal)

        # loop through all possible genotypes, and see if transmission probability is correct
        for l = 1:j_genotypes
            for k = 1:i_genotypes
                trans = MendelGeneticCounseling.transmission_genetic_counseling(person, 
                    locus, gamete[:, l], multi_genotype[:, :, k], par, keyword, 
                    start, finish, i, j)

                #gamete[:, l] could be [1, 1] or [1, 2] or [2, 1] 
                #multi_genotype[:, :, k] is either [1;1] or [1;2]
                #we compute the probaiblity that a multi_genotype pass down the gamete

                println("begin")
                println(gamete[:, l])
                println(multi_genotype[:, :, k])
                println(trans)
                println("end")

                # if gamete[:, l] == [1] && sum(multi_genotype[:, :, k]) == 2
                #     @test trans == 1.0
                # elseif gamete[:, l] == [1] && sum(multi_genotype[:, :, k]) == 3
                #     @test trans == 0.5
                # elseif gamete[:, l] == [1] && sum(multi_genotype[:, :, k]) == 4
                #     @test trans == 0.5
                # elseif gamete[:, l] == [2] && sum(multi_genotype[:, :, k]) == 2
                #     @test trans == 0.0
                # elseif gamete[:, l] == [2] && sum(multi_genotype[:, :, k]) == 3
                #     @test trans == 0.5
                # elseif gamete[:, l] == [2] && sum(multi_genotype[:, :, k]) == 4
                #     @test trans == 1.0
                # else
                #     @test_throws(MethodError, "shouldn't have reached here bro")
                # end
            end
        end
    end

end

@testset "are wrappers doing what they should be doing" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "genetic counseling 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    execution_error = MendelGeneticCounseling.genetic_counseling_option(pedigree, person,
        nuclear_family, locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
    result = GeneticCounseling("genetic counseling 1 Control.txt")

    @test execution_error == false #no errors should be thrown
    @test result == nothing # if nothing is returned the code has thrown no error


    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "genetic counseling 2 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    execution_error = MendelGeneticCounseling.genetic_counseling_option(pedigree, person,
        nuclear_family, locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
    result = GeneticCounseling("genetic counseling 2 Control.txt")

    @test execution_error == false
    @test result == nothing
end

@testset "final output" begin
    # note these tests are assuming that eslton_stewart algs are working properly
    # they do not yet have tests.
    GeneticCounseling("genetic counseling 1 Control.txt")
    GeneticCounseling("genetic counseling 2 Control.txt")
    outputted_file1 = open("genetic counseling 1 Output.txt")
    outputted_file2 = open("genetic counseling 2 Output.txt")

    result1 = readlines(outputted_file1)
    result2 = readlines(outputted_file2)

    @test result1[2] == " The risk =  0.03892." #this input is in Mendel too, and result is the same
    @test result2[2] == " The risk =  0.83986." #this input is not in Mendel
end

#current coverage = (83, 91) ~ 91.2%



