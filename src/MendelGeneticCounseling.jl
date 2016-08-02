"""
This module orchestrates a genetic counseling analysis.
"""
module MendelGeneticCounseling
#
# Required OpenMendel packages and modules.
#
using MendelBase
# using DataStructures                  # Now in MendelBase.
# using ModelConstruction               # Now in MendelBase.
# using ElstonStewartPreparation        # Now in MendelBase.
# using ElstonStewartEvaluation         # Now in MendelBase.
using Search
using SearchSetup
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
  keyword = set_keyword_defaults!(Dict{ASCIIString, Any}())
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
  keyword::Dict{ASCIIString, Any})
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
  loglikelihood = elston_stewart_loglikelihood(pedigree, person, locus, 
    parameter, instruction, keyword)
  risk = exp(pedigree.loglikelihood[1, 2] - pedigree.loglikelihood[2, 2])
  if risk > 1.0
    return execution_error = true
  end
  risk_string = @sprintf("%8.5f", risk)
  println(" The risk = $risk_string.")
  return execution_error = false
end # function genetic_counseling_option

end # module MendelGeneticCounseling

