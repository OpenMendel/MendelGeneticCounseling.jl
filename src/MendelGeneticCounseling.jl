__precompile__()

"""
This module orchestrates a genetic counseling analysis.
"""
module MendelGeneticCounseling
#
# Required OpenMendel packages and modules.
#
using MendelBase
using MendelSearch
#
# Required external modules.
#
using DataFrames
using Distributions
using LinearAlgebra
using Printf

export GeneticCounseling

"""
This is the wrapper function for the Genetic Counseling analysis option.
"""
function GeneticCounseling(control_file = ""; args...)
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println(" Genetic Counseling Analysis Option")
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "geneticcounseling")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "GeneticCounseling"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, person_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Insert mutation information into the locus data structure.
  # At most one type of mutation event is allowed.
  #
  process_mutation_keyword(locus, keyword)
  #
  # Check if SNP data were read.
  #
  if snpdata.snps != 0
    println(" \n\nERROR: This analysis does not use data from SNP files!\n")
  else
  #
  # Execute the specified analysis.
  #
    println(" \nAnalyzing the data.\n")
    execution_error = genetic_counseling_option(pedigree, person,
      nuclear_family, locus, locus_frame, phenotype_frame, person_frame,
      keyword)
    if execution_error
      println(" \n \nERROR: Mendel terminated prematurely!\n")
    else
      println(" \n \nMendel's analysis is finished.\n")
    end
  end
  #
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing
end # function GeneticCounseling

"""
This option requires a numerator pedigree (number 1) and a
denominator pedigree (number 2). It computes a risk equal to
the ratio of the likelihoods of these two pedigrees.
The function returns an execution error indicator.
"""
function genetic_counseling_option(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus, locus_frame::DataFrame,
  phenotype_frame::DataFrame, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any})
  #
  # Unless there is mutation in a counseling problem,
  # eliminate genotypes and lump alleles.
  #
  keyword["eliminate_genotypes"] = true
  keyword["lump_alleles"] = true
  #
  # Define the parameter data structure.
  #
  parameter = set_parameter_defaults(keyword)
  #
  # Fetch the instructions for conducting the Elston-Stewart algorithm.
  #
  (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
    person, nuclear_family, locus, keyword)
  if any(elston_stewart_count .>  keyword["complexity_threshold"])
    println("One or more pedigrees exceeded the complexity threshold.")
    println(io, "One or more pedigrees exceeded the complexity threshold.")
    return execution_error = true
  end
  #
  # Evaluate the risk via the Elston-Stewart algorithm.
  #
  loglikelihood = elston_stewart_loglikelihood(penetrance_genetic_counseling,
    prior_genetic_counseling, transmission_genetic_counseling,
    pedigree, person, locus, parameter, instruction, person_frame, keyword)
  risk = exp(pedigree.loglikelihood[1, 2] - pedigree.loglikelihood[2, 2])
  if risk > 1.0
    return execution_error = true
  end
  #
  # Output the result, first to the screen, then to a file.
  #
  risk_string = @sprintf("%8.5f", risk)
  println(" The risk = $risk_string.")
  if keyword["output_file"] != ""
    output_unit = keyword["output_unit"]
    println(output_unit, "\n The risk = $risk_string.\n \n")
  end
  return execution_error = false
end # function genetic_counseling_option

"""
Supply a penetrance for individual i.
"""
function penetrance_genetic_counseling(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64}, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any}, startlocus::Int, finishlocus::Int, i::Int)

  pen = 1.0
  for l = startlocus:finishlocus
    allele1 = multi_genotype[1, l]
    allele2 = multi_genotype[2, l]
    p = 1.0 # for reduced penetrance let p depend on loc
    loc = locus.model_locus[l]
    if loc == locus.mutable_locus 
#       x = person_frame[i, :LogCPK]
#       if !ismissing(x)
#         glm_scale = 0.3255
#         glm_mean = 1.835
#         if allele1 == 2 || allele2 == 2
#           glm_mean = glm_mean + 0.265  
#         else
#           glm_mean = glm_mean - 0.265 
#         end
#         dist = Normal(glm_mean, glm_scale)
#         p = pdf.(dist, x)
#       end
    end
    pen = p * pen
  end
  return pen
end # function penetrance_genetic_counseling

"""
Supply a prior probability for founder i.
"""
function prior_genetic_counseling(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64}, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any}, startlocus::Int, finishlocus::Int, i::Int)

  prior_prob = 1.0
  for l = startlocus:finishlocus
    loc = locus.model_locus[l]
    allele = multi_genotype[1, l]
    frequency = dot(person.admixture[i, :], locus.frequency[loc][:, allele])
    prior_prob = prior_prob * frequency
    if !locus.xlinked[loc] || !person.male[i]
      allele = multi_genotype[2, l]
      frequency = dot(person.admixture[i, :], locus.frequency[loc][:, allele])
      prior_prob = prior_prob * frequency
    end
  end
  return prior_prob
end # function prior_genetic_counseling

"""
Supply the transmission probability that a parent i with a particular
genotype transmits a particular gamete to his or her child j.
This version of the transmission function allows mutation.
"""
function transmission_genetic_counseling(person::Person, locus::Locus,
  gamete::Vector{Int}, multi_genotype::Matrix{Int}, par::Vector{Float64},
  person_frame::DataFrame, keyword::Dict{AbstractString, Any},
  startlocus::Int, finishlocus::Int, i::Int, j::Int)
  #
  # Store the sex of the parent as an integer; 1 for female and 2 for male.
  #
  i_sex = person.male[i] + 1 
  xlinked = locus.xlinked[locus.model_locus[startlocus]]
  rate = locus.mutation_rate[i_sex]
  #
  # For male to male inheritance at an x-linked locus,
  # set the transmission probability equal to 1.
  # For male to female transmission, account for matching and mutation.
  #
  trans = 1.0
  if xlinked && person.male[i]
    if person.male[j]
      return trans
    else
      for l = startlocus:finishlocus
        (k, m) = (multi_genotype[2, l], gamete[l])
        if locus.model_locus[l] !== locus.mutable_locus
          if m != k
            return trans = 0.0
          end
        else
          if k == locus.donor_allele
            if m == locus.recipient_allele
              trans = trans * rate
            elseif m == k
              trans = trans * (1.0 - rate)
            else
              return trans = 0.0
            end
          end
        end
      end
    end
    return trans
  end
  #
  # Consider autosomal transmission or maternal X-linked transmission.
  #
  both = trues(b)
  source = ones(Int, b)
  for l = startlocus:finishlocus
    (k, m1, m2) = (gamete[l], multi_genotype[1, l], multi_genotype[2, l])
    match1 = k == m1
    match2 = k == m2
    #
    # Allow for a mutational match.
    #
    if locus.model_locus[l] == locus.mutable_locus
      if k == locus.recipient_allele
        match1 = match1 || m1 == locus.donor_allele 
        match2 = match2 || m2 == locus.donor_allele 
      end
    end
    #
    # Check whether either the first or second parental gene at the 
    # current locus matches or is mutable to the gamete gene at this 
    # locus. If not, then return with 0 for the transmission probability.
    #
    if !match1 && !match2
      return trans = 0.0
    end
    #
    # Initialize the source vector.
    #
    if !match1
      source[l] = 2
    end
    both[l] = match1 && match2
  end
  #
  # Generate the next gamete from the parental genes which is identical
  # or mutable to the child's gamete.  These are taken in reverse dictionary
  # order as strings from the two letter alphabet with letters 1 and 2.
  #
  (trans, found, finished, first_gamete) = (0.0, false, false, true)
  while !finished
    if !first_gamete
      for l = startlocus:finishlocus
        found = false
        if both[l] 
          if source[l] == 1 
            source[l] = 2
            found = true
            break
          else
            source[l] = 1
          end
        end
      end
      if !found
        finished = true 
        break
      end
    end
    first_gamete = false
    if !finished
      #
      # Initialize the transmission probability for this parental gamete.
      #
      if startlocus == 1 
        r = 0.5
      else
        r = 1.0
      end
      #
      # Compute the remaining factors of the transmission probability.
      # Start by including mutation.
      #
      for l = startlocus:finishlocus
        if locus.model_locus[l] == locus.mutable_locus
          if multi_genotype[source[l], l] == locus.donor_allele 
            if gamete[l] == locus.recipient_allele 
              r = r * rate
            else
              r = r * (1.0 - rate)
            end
          end
        end
        #
        # Pick up the recombination fractions.
        #
        if l < finishlocus 
          if source[l] == source[l + 1] 
            r = r * (1.0 - locus.theta[i_sex, l])
          else
            r = r * locus.theta[i_sex, l]
          end
        end
      end
      trans = trans + r
    end
  end
  return trans
end # function transmission_genetic_counseling
#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

"""
Checks and interprets the mutation keyword.
"""
function process_mutation_keyword(locus::Locus,
  keyword::Dict{AbstractString, Any})
#
  println("mutation info = ",keyword["mutation"])
  value = split(keyword["mutation"], ',', keepempty = false)
  if length(value) == 0
    locus.mutable_locus = 0
    return
  elseif length(value) != 5
    throw(ArgumentError(
        "The keyword mutation must take exactly five values.\n \n"))
  end
  loc = something(findfirst(isequal(value[1]), locus.name), 0)
  if loc == 0
    throw(ArgumentError(
        "Not a valid name for a mutating locus.\n \n"))
  else    
    locus.mutable_locus = loc
  end
  i = something(findfirst(isequal(value[2]), locus.allele_name[loc]), 0) 
  if i == 0
    throw(ArgumentError(
        "Not a valid name for a mutating allele.\n \n"))
  else  
    locus.donor_allele = i
  end
  i = something(findfirst(isequal(value[3]), locus.allele_name[loc]), 0) 
  if i == 0
    throw(ArgumentError(
        "Not a valid name for an allele mutated to.\n \n"))
  else
    locus.recipient_allele = i  
  end
  v4 = parse(Float64, value[4])
  if v4 < 0.0 || v4 > 1.0
     throw(ArgumentError(
        "Not a valid female mutation rate.\n \n"))
  else
    locus.mutation_rate[1] = v4 
  end
  v5 = parse(Float64, value[5])
  if v5 < 0.0 || v5 > 1.0
     throw(ArgumentError(
        "Not a valid fmale mutation rate.\n \n"))
  else
    locus.mutation_rate[2] = v5 
  end
end

end # module MendelGeneticCounseling
