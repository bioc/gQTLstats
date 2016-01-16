mixedVCFtoSnpMatrix = function(vcf, preferGT=TRUE) {
#
# current behavior of genotypeToSnpMatrix is to use either GT
# or one of GL, GP, PL to obtain bytecodes for call or allelic dose
# the "or" is exclusive
# 
# this function will obtain the GT-based and probability-based
# results, substituting the latter for the former according to the
# setting of preferGT
#
  r1 = genotypeToSnpMatrix(vcf, uncertain=FALSE)
  r2 = genotypeToSnpMatrix(vcf, uncertain=TRUE)
  NACODE = as.raw(0)
  matOfCalls = r1$genotypes@.Data
  matOfProbs = r2$genotypes@.Data
# substitute codes from probabilities into missing calls
  if (preferGT) { # only use probability-based codes for missing calls
    matOfCalls[ matOfCalls == NACODE ] = matOfProbs[ matOfCalls == NACODE ]
    }
  else {  # use all nonmissing prob-based codes
    matOfCalls[ matOfProbs != NACODE ] = matOfProbs[ matOfProbs != NACODE ]
    }
  r1$genotypes@.Data = matOfCalls
  r1
}
