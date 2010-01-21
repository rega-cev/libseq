ADD_EXECUTABLE(aafastaread src/AAFastaRead.C)
ADD_EXECUTABLE(ntfastaread src/NTFastaRead.C)
ADD_EXECUTABLE(tajimad src/TajimaD.C)
ADD_EXECUTABLE(codonalign src/CodonAlign.C)
ADD_EXECUTABLE(nmwtest src/NeedlemanWunshTest.C)
ADD_EXECUTABLE(genetic_diversity src/GeneticDiversity.C)
ADD_EXECUTABLE(subscore src/SubstitutionScoreRegression.C)

TARGET_LINK_LIBRARIES(subscore seq)
TARGET_LINK_LIBRARIES(aafastaread seq)
TARGET_LINK_LIBRARIES(ntfastaread seq)
TARGET_LINK_LIBRARIES(nmwtest seq)
TARGET_LINK_LIBRARIES(tajimad seq)
TARGET_LINK_LIBRARIES(codonalign seq)
TARGET_LINK_LIBRARIES(genetic_diversity seq)
INCLUDE_DIRECTORIES(
  ${SEQ_SOURCE_DIR}/src/sequence
	${SEQ_SOURCE_DIR}/src/evolution
	${SEQ_SOURCE_DIR}/src/algorithm
)
