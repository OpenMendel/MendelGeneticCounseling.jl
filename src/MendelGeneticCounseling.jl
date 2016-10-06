"""
This module orchestrates a genetic counseling analysis.
"""
module MendelGeneticCounseling
#
# Required OpenMendel packages and modules.
#
using MendelBase                        # From package MendelBase.
# using DataStructures                  # Now in MendelBase.
# using ModelConstruction               # Now in MendelBase.
# using ElstonStewartPreparation        # Now in MendelBase.
# using ElstonStewartEvaluation         # Now in MendelBase.
using SearchSetup                       # From package Search.
#
# Required external modules.
#
using DataFrames                        # From package DataFrames.

export GeneticCounseling

"""
This is the wrapper function for the Genetic Counseling analysis option.
"""
function GeneticCounseling(control_file = ""; args...)

  const GENETIC_COUNSELING_VERSION :: VersionNumber = v"0.1.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println(" Genetic Counseling analysis option")
  println("        version ", GENETIC_COUNSELING_VERSION )
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
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Execute the specified analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = genetic_counseling_option(pedigree, person,
    nuclear_family, locus, locus_frame, phenotype_frame, pedigree_frame,
    keyword)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
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
  phenotype_frame::DataFrame, pedigree_frame::DataFrame,
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
    pedigree, person, locus, parameter, instruction, keyword)
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
  multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int)

  pen = 1.0
  for l = start:finish
    allele1 = multi_genotype[1, l]
    allele2 = multi_genotype[2, l]
    loc = locus.model_locus[l]
    p = 1.0 # for reduced penetrance let p depend on loc
    pen = p * pen
  end
  return pen
end # function penetrance_genetic_counseling

"""
Supply a prior probability for founder i.
"""
function prior_genetic_counseling(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int)

  prior_prob = 1.0
  for l = start:finish
    loc = locus.model_locus[l]
    allele = multi_genotype[1, l]
    frequency = dot(vec(person.admixture[i, :]),
                    vec(locus.frequency[loc][:, allele]))
    prior_prob = prior_prob * frequency
    if !locus.xlinked[loc] || !person.male[i]
      allele = multi_genotype[2, l]
      frequency = dot(vec(person.admixture[i, :]),
                      vec(locus.frequency[loc][:, allele]))
      prior_prob = prior_prob * frequency
    end
  end
  return prior_prob
end # function prior_genetic_counseling

"""
Supply the transmission probability that a parent i with a particular
genotype transmits a particular gamete to his or her child j.
"""
function transmission_genetic_counseling(person::Person, locus::Locus,
  gamete::Vector{Int}, multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int, j::Int)
  #
  # For male to male inheritance at an x-linked locus,
  # set the transmission probability equal to 1.
  #
  loc = locus.model_locus[start]
  xlinked = locus.xlinked[loc]
  if xlinked && person.male[i] && person.male[j]
    return 1.0
  end
  #
  # Allow for linked markers but no mutation. Start by storing an indicator
  # of the sex of the parent.
  #
  if person.male[i]
    i_sex = 2
  else
    i_sex = 1
  end
  #
  # Reduce the computations by considering only the heterozygous loci.
  # Use Trow's formula to express the recombination fraction
  # between two heterozygous loci in terms of the recombination
  # fractions between the adjacent loci that separate them.
  # Set the logical variable found to true when the first heterozygous
  # parental locus is found. Phase records the phase of the most
  # recent heterozygous parental locus.
  #
  trans = 1.0
  found = false
  phase = true
  r = 0.5
  for l = start:finish
    match1 = multi_genotype[1, l] == gamete[l]
    match2 = multi_genotype[2, l] == gamete[l]
    #
    # Check whether either the first or second parental gene at
    # the current locus matches the gamete gene at this locus.
    # If not, then return with 0 for the transmission probability.
    #
    if !match1 && !match2
      return 0.0
    end
    #
    # Check whether the current locus is heterozygous.
    #
    if match1 != match2
      if found
        if phase == match1
          trans = trans * (0.5 + r)
        else
          trans = trans * (0.5 - r)
        end
      else
        found = true
        if start == 1 || start == finish
          trans = 0.5
        else
          trans = 1.0
        end
      end
      phase = match1
      r = 0.5
    end
    if found && l < finish
      r = r * (1.0 - 2.0 * locus.theta[i_sex, l])
    end
  end
  if !found; trans = 1.0; end
  return trans
end # function transmission_genetic_counseling

end # module MendelGeneticCounseling

