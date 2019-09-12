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
using LinearAlgebra
using Printf
using Distributions
using CSV
###
## adding in the inverse links so I can use the appropriate distribution
###
  include("apply_inverse_link_new.jl")
  export LogLink, IdentityLink, SqrtLink, ProbitLink, LogitLink, InverseLink, CauchitLink, CloglogLink
######
### try adding in distributions so I can bring in the appropriate data
  include("apply_dist.jl")
  export BinomialDist, NegativeBinomialDist, PoissonDist, ExponentialDist, GammaDist, InverseGaussianDist, LogisticDist, LognormalDist, NormalDist
######
export GeneticCounseling

"""
This is the wrapper function for the Genetic Counseling analysis option.

"""
######
# The first task is to find the variables to make the mean vector:
find_variables(x) = find_variables!(Symbol[], x) #this is so that we can call the function without the exclamation
function find_variables!(var_names, x::Number) #if the variable name is a number then just return it without doing anything
    return var_names
end
function find_variables!(var_names, x::Symbol) # if the variable name is a symbol then push it to the list of var_names because its a name of a column
    push!(var_names, x)
end

function find_variables!(var_names, x::Expr) # if the variable is a expression object then we have to crawl through each argument of the expression
    # safety checking
    if x.head == :call  # check for + symbol bc we are summing linear combinations within each expression argument
      # pass the remaining expression
      for argument in x.args[2:end] # since the first argument is the :+ call
        find_variables!(var_names, argument) #recursively find the names given in each argument so check if number, if symbol, if expression etc.again agian
    end
end
return var_names
end

######
#######
function search_variables!(x::Expr, var::Symbol)
    for i in eachindex(x.args)
        if x.args[i] == var # if the argument is one of the variables given then just put it in the right format df[:x1]
###            x.args[i] = Meta.parse(string(:input_data_from_user,"[", ":", var, "]"))
            x.args[i] = Meta.parse(string(:input_data_from_user,"[:, :", var, "]"))
        elseif x.args[i] isa Expr # else if the argument is an expression (i.e not a varaible (symbol) or a number) then
            search_variables!(x.args[i], var) #go through this function recursively on each of the arguments of the expression object
        end
    end
    return x
end

function search_variables!(x::Expr, vars...) # this is for when you have more than one variable name found in the string
    for var in vars #goes through each of the variables in the vector vars
        x = search_variables!(x, var) #runs the recursion on each variable in vars
    end
    return x
end
#######
##
function mean_formula(user_formula_string::String, df::DataFrame)
    global input_data_from_user = df

    users_formula_expression = Meta.parse(user_formula_string)
    if(users_formula_expression isa Expr)
        found_markers = find_variables(users_formula_expression) #store the vector of symbols of the found variables

        dotted_args = map(Base.Broadcast.__dot__, users_formula_expression.args) # adds dots to the arguments in the expression
        dotted_expression = Expr(:., dotted_args[1], Expr(:tuple, dotted_args[2:end]...)) #reformats the exprssion arguments by changing the variable names to tuples of the variable names to keep the dot structure of julia

        julia_interpretable_expression = search_variables!(dotted_expression, found_markers...) #gives me the julia interpretable exprsesion with the dataframe provided

        mean_vector = eval(Meta.parse(string(julia_interpretable_expression))) #evaluates the julia interpretable expression on the dataframe provided
    else
        mean_vector = [users_formula_expression for i in 1:size(df, 1)]
    end
    return mean_vector
end
######
####### 
#######
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
  ####
  # keywords needed for Genetic Counseling are glm_mean, inverse_link
  # distribution
    keyword["glm_mean"] = "0.0"
    keyword["glm_link"] = "IdentityLink"
    keyword["glm_response"] = "NormalDist"
    keyword["glm_trait"] = "Affected"
    keyword["glm_trials"]=1
    keyword["glm_scale"]=1.0
    keyword["penetrance_file"] = ""
    keyword["penetrance"] = ""
    keyword["riskfactor"]= ""
    keyword["disease_status"]=""
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
  # Check if SNP data were read.
  #
  if snpdata.snps != 0
    println(" \n\nERROR: This analysis does not use data from SNP files!\n")
  else
  #
  # Execute the specified analysis.
  #
####
## add in two columns for the alleles of an individual at a putative locus
#####
##
### bringing in a penetrance dictionary
###
#
#Starting to make penetrance function conditional of existing of penetrance file
#
   penetrancefile = string(keyword["penetrance_file"])
#
     if keyword["penetrance_file"] == ""
        println(" no penetrance file")
     else
        pen_df = CSV.File(penetrancefile)  |> DataFrame
        names_pen = names(pen_df)
        names_ped = names(pedigree_frame)
        names_risk = intersect(names_pen, names_ped)
        risk_factors = length(names_risk)
#        status = string(keyword["disease_status"])
        println("this problem has $risk_factors factors called $names_risk")
#
# Make penetrance dataframe into a dictionary
#
        pen_dict = Dict{Tuple, Array{Float64, 1}}()
        for j = 1:size(pen_df, 1)
           pen_dict[Tuple(pen_df[j, 4:end])] = convert(Array{Float64,1},pen_df[j,1:3])
        end
        keyword["penetrance"]=pen_dict
        keyword["riskfactor"]=names_risk
     end

        q = size(pedigree_frame,1)
##    pedigree_frame.allele1 =  missings(Int,q)
##    pedigree_frame.allele2 = missings(Int,q)
    pedigree_frame[!, :allele1] = missings(Int,q)
    pedigree_frame[!, :allele2] = missings(Int,q)

    println(" \nAnalyzing the data.\n")
    execution_error = genetic_counseling_option(pedigree, person,
      nuclear_family, locus, locus_frame, phenotype_frame, pedigree_frame,
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
    pedigree, person, locus, parameter, instruction, pedigree_frame, keyword)
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
  multi_genotype::Matrix{Int}, par::Vector{Float64}, pedigree_frame::DataFrame,
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int)
#####
### bringing in the distribution from a keyword
### need to figure out how to handle converting from a string into
### a distribution, mean and possibly a scale
####
#   glmresponse() = Normal()
  pen = 1.0
  glmtrial = keyword["glm_trials"]
  glmscale = keyword["glm_scale"]
# pick my trait 
#
  glmtrait = Symbol(keyword["glm_trait"])
  status = Symbol(keyword["disease_status"])
#
#
# Bring in the Penetrance information
#
   penetrancefile = string(keyword["penetrance_file"])
#
# Checking on when penetrance file name changes so can write conditional
#
    if keyword["penetrance_file"] != ""
       disease_entry = pedigree_frame[i,status]
       identifier = pedigree_frame[i,:Person]
#       println(" disease status is for individual $identifier is $disease_entry")
    end
#
# Calculate penetrance for each individual at each valid genotype
#
  for l = start:finish
    allele1 = multi_genotype[1, l]
    allele2 = multi_genotype[2, l]
    loc = locus.model_locus[l]
##    pedigree_frame.allele1[i]=allele1
##    pedigree_frame.allele2[i]=allele2
    pedigree_frame[i, :allele1]=allele1
    pedigree_frame[i, :allele2]=allele2
    single_frame = DataFrame(pedigree_frame[i, :])
    if keyword["penetrance_file"] ==""
         glmmean = string(keyword["glm_mean"])
         glminvlink = Symbol(keyword["glm_link"])
         w = mean_formula(glmmean,single_frame)   
         finv_link = getfield(Main,glminvlink)
         transmu = apply_inverse_link(w,finv_link())
         glmresponse = Symbol(keyword["glm_response"])
         fmydist = getfield(Main,glmresponse)
#  this version allows the user to specify the name of the trait to use
           if ismissing(single_frame[1,glmtrait]) 
#          if single_frame[1,glmtrait]<0.0
            p = 1.0
         else
            penout = pdf.(apply_dist(glmtrial,glmscale,transmu,fmydist()),single_frame[1,glmtrait])
            p = penout[1]
         end
#    println("pdf of the penetrance is $p")
     pen = p * pen
#    println("Individual $identifier has penetrance $pen when their genotype is $allele1 - $allele2")
    else
        risknames = keyword["riskfactor"]
          icount = allele1+allele2-1
#          myentry=pen_dict[(pedigree_frame[i,Symbol(names_risk[1])],pedigree_frame[i,Symbol(names_risk[2])])][icount]
           myentry=keyword["penetrance"][(pedigree_frame[i,Symbol(risknames[1])],pedigree_frame[i,Symbol(risknames[2])])][icount]
# keyword["penetrance"]
          if disease_entry == 1 
              pen = pen * myentry
          elseif disease_entry == 0
              pen = pen * (1.0 - myentry)
          else
              pen = 1.0      
          end 
    end
#      println("Individual $identifier has penetrance $pen when their genotype is $allele1 - $allele2 and their disease status is $disease_entry")
  end
  return pen
  println(" check value $pen")
end # function penetrance_genetic_counseling
"""
Supply a prior probability for founder i.
"""
function prior_genetic_counseling(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64}, pedigree_frame::DataFrame,
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int)

  prior_prob = 1.0
  for l = start:finish
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
"""
function transmission_genetic_counseling(person::Person, locus::Locus,
  gamete::Vector{Int}, multi_genotype::Matrix{Int}, par::Vector{Float64},
  pedigree_frame::DataFrame, keyword::Dict{AbstractString, Any}, start::Int, 
  finish::Int, i::Int, j::Int)
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
#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module MendelGeneticCounseling
