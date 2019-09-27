package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ahrd.controller.AHRD;
import ahrd.model.BlastResult;
import ahrd.model.CompetitorAnnotation;
import ahrd.model.Protein;

public class EvaluatorOutputWriter extends TsvOutputWriter {

	protected BufferedWriter hrdScoresWriter;
	private String columnNamesLine = "";
	private Map<Protein, String> proteinOutputLines = new HashMap<Protein, String>();
	private String averagesLine = "";
	private String coveragesLine = "";
	private final String separator = "\t";

	public EvaluatorOutputWriter(Collection<Protein> proteins) {
		super(proteins);
	}

	public void writeOutput() throws IOException {
		/**
		 * Prepare the output
		 */
		prepareOutput();
		
		/**
		 * Write the output to file
		 */
		BufferedWriter bw = new BufferedWriter(new FileWriter(getSettings().getPathToOutput()));
		// Header
		bw.write("# AHRD-Version " + AHRD.VERSION + " (Evaluator)\n");
		bw.write("\n");
		// Column names:
		bw.write(columnNamesLine + "\n");
		// Protein data
		for (Protein prot : getProteins()) {
			bw.write(proteinOutputLines.get(prot) + "\n");
		}
		// Summary
		if (getSettings().doWriteEvaluationSummary()) {
			bw.write("\n");
			bw.write(averagesLine + "\n");
			bw.write(coveragesLine);
		}
		// Cleanup
		bw.close();

		/**
		 *  If AHRD is requested to write out the AHRD-Score of each BlastHit's Description, do so into another file
		 */
		if (getSettings().doWriteHRDScoresToOutput()) {
			writeHRDScoresOutputHeader();
			for (Protein prot : getProteins()) {
				writeHrdScoresOutput(prot);
			}
			this.hrdScoresWriter.close();
		}
	}
	
	private void prepareOutput() {
		// AHRD's most basic columns
		columnNamesLine += "Protein-Accession" + separator + "Blast-Hit-Accession" + separator + "AHRD-Quality-Code" + separator + "Human-Readable-Description";
		for (Protein prot : getProteins()) {
			proteinOutputLines.put(prot, buildDescription(prot));
		}
		averagesLine += "Average" + separator + separator + separator;
		coveragesLine += "Coverage" + separator + separator + separator;

		// AHRD's GO terms:
		if (getSettings().hasGeneOntologyAnnotations()) {
			columnNamesLine += separator + "Gene-Ontology-Terms";
			for (Protein prot : getProteins()) {
				proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + combineGoTermStrings(prot.getGoResults()));
			}
			averagesLine += separator;
			coveragesLine += separator;
		}
		
		// AHRD's description evaluation
		columnNamesLine += separator + "HRD-Length" + separator + "Ground-Truth-Description" + separator + "Ground-Truth-Description-Length" + separator + "Evaluation-Score";
		int sumDescriptionLength = 0;
		int covDescriptLength = 0;
		int sumGroundTruthDescriptionLength = 0;
		int covGroundTruthDescriptionLength = 0;
		double sumEvaluationScore = 0.0;
		int covEvaluationScore = 0;
		for (Protein prot : getProteins()) {
			String descriptionEvaluationLine = "";
			if (prot.getEvaluationScoreCalculator().getEvalutionScore() != null) {
				// HRD-Length
				if (prot.getDescriptionScoreCalculator().getHighestScoringBlastResult() != null) {
					int descriptionLength = prot.getDescriptionScoreCalculator().getHighestScoringBlastResult().getEvaluationTokens().size();
					descriptionEvaluationLine += descriptionLength;
					sumDescriptionLength += descriptionLength;
					covDescriptLength++;
				}
				else {
					descriptionEvaluationLine += "NaN";
				}
				// Ground truth description
				descriptionEvaluationLine += separator + prot.getEvaluationScoreCalculator().getGroundTruthDescription().getDescription();
				// Ground truth description length
				int groundTruthDescriptionLength = prot.getEvaluationScoreCalculator().getGroundTruthDescription().getTokens().size();
				if (groundTruthDescriptionLength > 0) {
					descriptionEvaluationLine += separator + groundTruthDescriptionLength;
					sumGroundTruthDescriptionLength += groundTruthDescriptionLength;
					covGroundTruthDescriptionLength++;
				} else {
					descriptionEvaluationLine += separator + "NaN";
				}
				// AHRD's description score
				Double evaluationScore = prot.getEvaluationScoreCalculator().getEvalutionScore().getScore();
				descriptionEvaluationLine += separator + FRMT.format(evaluationScore);
				if (!evaluationScore.isNaN()) {
					sumEvaluationScore += evaluationScore;
					covEvaluationScore++;
				}
			} else {
				descriptionEvaluationLine += separator + separator + separator;
			}
			proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + descriptionEvaluationLine);
		}
		averagesLine += separator + FRMT.format((double)sumDescriptionLength/covDescriptLength) + separator + separator + FRMT.format((double)sumGroundTruthDescriptionLength/covGroundTruthDescriptionLength) + separator + FRMT.format(sumEvaluationScore/covEvaluationScore);
		coveragesLine += separator + FRMT.format((double)covDescriptLength/getProteins().size()) + separator + separator + FRMT.format((double)covGroundTruthDescriptionLength/getProteins().size()) + separator + FRMT.format((double)covEvaluationScore/getProteins().size());
		
		// AHRD's description evaluation details
		if (getSettings().doWriteFscoreDetailsToOutput()) {
			columnNamesLine += separator +  "Diff-to-bestCompetitor" + separator + "Precision" + separator + "Recall";
			double sumEvalScoreMinBestCompScore = 0.0;
			int covEvalScoreMinBestCompScore = 0;
			double sumEvaluationScorePrecision = 0.0;
			int covEvaluationScorePrecision = 0;
			double sumEvaluationScoreRecall = 0.0;
			int covEvaluationScoreRecall = 0;
			for (Protein prot : getProteins()) {
				String descriptionEvaluationDetailsLine = "";
				if (prot.getEvaluationScoreCalculator().getEvalutionScore() != null) {
					Double evalScoreMinBestCompScore = prot.getEvaluationScoreCalculator().getEvalScoreMinBestCompScore();
					descriptionEvaluationDetailsLine += FRMT.format(evalScoreMinBestCompScore);
					if (!evalScoreMinBestCompScore.isNaN()) {
						sumEvalScoreMinBestCompScore += evalScoreMinBestCompScore;
						covEvalScoreMinBestCompScore++;
					}
					Double evaluationScorePrecision = prot.getEvaluationScoreCalculator().getEvalutionScore().getPrecision();
					descriptionEvaluationDetailsLine += separator + FRMT.format(evaluationScorePrecision);
					if (!evaluationScorePrecision.isNaN()) {
						sumEvaluationScorePrecision += evaluationScorePrecision;
						covEvaluationScorePrecision++;
					}
					Double evaluationScoreRecall = prot.getEvaluationScoreCalculator().getEvalutionScore().getRecall();
					descriptionEvaluationDetailsLine += separator + FRMT.format(evaluationScoreRecall);
					if (!evaluationScoreRecall.isNaN()) {
						sumEvaluationScoreRecall += evaluationScoreRecall;
						covEvaluationScoreRecall++;
					}
				} else {
					descriptionEvaluationDetailsLine += separator + separator + separator;
				}
				proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + descriptionEvaluationDetailsLine);
			}
			averagesLine += separator + sumEvalScoreMinBestCompScore/covEvalScoreMinBestCompScore + separator + sumEvaluationScorePrecision/covEvaluationScorePrecision + separator + sumEvaluationScoreRecall/covEvaluationScoreRecall;
			coveragesLine += separator + (double)covEvalScoreMinBestCompScore/getProteins().size() + separator + (double)covEvaluationScorePrecision/getProteins().size() + separator + (double)covEvaluationScoreRecall/getProteins().size();
		}
		
		if (getSettings().getWriteDescriptionSubScoresToOutput()) {
			columnNamesLine += separator + "Sum(Token-Scores)" + separator + "TokenHighScore" + separator + "Correction-Factor" + separator + "Lexical-Score" + separator + "RelativeBlastScore" + separator + "DescriptionScore";
			for (Protein prot : getProteins()) {
				BlastResult highestScoringBlastResult = prot.getDescriptionScoreCalculator().getHighestScoringBlastResult();
				// Found a high scoring description?
				if (highestScoringBlastResult == null) {
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + "NaN" + separator + "NaN" + separator + "NaN" + separator + "NaN" + separator + "NaN" + separator + "NaN");
				} else {
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(prot.getTokenScoreCalculator().sumOfAllTokenScores(highestScoringBlastResult))); 
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(prot.getTokenScoreCalculator().getTokenHighScore()));
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(prot.getLexicalScoreCalculator().correctionFactor(highestScoringBlastResult)));
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(prot.getLexicalScoreCalculator().lexicalScore(highestScoringBlastResult)));
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(prot.getDescriptionScoreCalculator().relativeBlastScore(highestScoringBlastResult)));
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestScoringBlastResult.getDescriptionScore()));
				}
			}
			averagesLine += separator + separator + separator + separator + separator + separator;
			coveragesLine += separator + separator + separator + separator + separator + separator;
		}
		
		// Best blast hits
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			for (String blastDb : getSettings().getBlastDatabases()) {
				if (blastDb != null && !blastDb.equals("")) {
					// Best blast hits accession and description
					columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-Accession-and-Description";
					for (Protein prot : getProteins()) {
						if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
							BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + "\"" + bestBr.getAccession() + " " + bestBr.getDescription() + "\"");
						} else {
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
						}
					}
					averagesLine += separator;
					coveragesLine += separator;
					if (getSettings().isInEvaluationMode()) {
						// Best blast hits description evaluation
						columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-Description-Length" + separator + "Best-BlastHit-against-" + blastDb + "-Description-Evaluation-Score";
						int sumBestBlastDescriptionLength = 0;
						int covBestBlastDescriptionLength = 0;
						Double sumBestBlastDescriptionScore = 0.0;
						int covBestBlastDescriptionScore = 0;
						for (Protein prot : getProteins()) {
							if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
								BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
								int bestBlastDescriptionLength = bestBr.getEvaluationTokens().size();
								if (bestBlastDescriptionLength > 0) {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + bestBlastDescriptionLength);
									sumBestBlastDescriptionLength += bestBlastDescriptionLength;
									covBestBlastDescriptionLength++;
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + "NaN");
								}
								Double bestBlastDescriptionScore = bestBr.getEvaluationScore().getScore();
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastDescriptionScore));
								if (!bestBlastDescriptionScore.isNaN()) {
									sumBestBlastDescriptionScore += bestBlastDescriptionScore;
									covBestBlastDescriptionScore++;
								}
							} else {
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
							}
						}
						averagesLine += separator + FRMT.format((double)sumBestBlastDescriptionLength/covBestBlastDescriptionLength) + separator + FRMT.format(sumBestBlastDescriptionScore/covBestBlastDescriptionScore);
						coveragesLine += separator + FRMT.format((double)covBestBlastDescriptionLength/getProteins().size()) + separator + FRMT.format((double)covBestBlastDescriptionScore/getProteins().size());
						if (getSettings().doWriteFscoreDetailsToOutput()) {
							columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-Description-Evaluation-Score-Precision\tBest-BlastHit-against-" + blastDb + "-Description-Evaluation-Score-Recall";
							double sumBestBlastDescriptionScorePrecision = 0.0;
							int covBestBlastDescriptionScorePrecision = 0;
							double sumBestBlastDescriptionScoreRecall = 0.0;
							int covBestBlastDescriptionScoreRecall = 0;
							for (Protein prot : getProteins()) {
								if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
									BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
									Double bestBlastDescriptionScorePrecision = bestBr.getEvaluationScore().getPrecision();
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastDescriptionScorePrecision));
									if (!bestBlastDescriptionScorePrecision.isNaN()) {
										sumBestBlastDescriptionScorePrecision += bestBlastDescriptionScorePrecision;
										covBestBlastDescriptionScorePrecision++;
									}
									Double bestBlastDescriptionScoreRecall = bestBr.getEvaluationScore().getRecall();
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastDescriptionScoreRecall));
									if (!bestBlastDescriptionScoreRecall.isNaN()) {
										sumBestBlastDescriptionScoreRecall += bestBlastDescriptionScoreRecall;
										covBestBlastDescriptionScoreRecall++;
									}
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
								}
							}
							averagesLine += separator + FRMT.format(sumBestBlastDescriptionScorePrecision/covBestBlastDescriptionScorePrecision) + separator + FRMT.format(sumBestBlastDescriptionScoreRecall/covBestBlastDescriptionScoreRecall);
							coveragesLine += separator + FRMT.format((double)covBestBlastDescriptionScorePrecision/getProteins().size()) + separator + FRMT.format((double)covBestBlastDescriptionScoreRecall/getProteins().size()); 
						}
						// Best blast hits go terms
						if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
							columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations";
							for (Protein prot : getProteins()) {
								if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
									BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + combineGoTermsToString(bestBr.getGoAnnotations()));
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
								}
							}
							averagesLine += separator;
							coveragesLine += separator;
							// Best blast hits go term evaluation
							if (getSettings().doCalculateSimpleGoF1Scores()) {
								columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-Simple-F-Score";
								double sumBestBlastGoAnnotationsSimpleScore = 0.0;
								int covBestBlastGoAnnotationsSimpleScore = 0;
								for (Protein prot : getProteins()) {
									if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
										BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
										Double bestBlastGoAnnotationSimpleScore = bestBr.getSimpleGoAnnotationScore().getScore();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationSimpleScore));
										if (!bestBlastGoAnnotationSimpleScore.isNaN()) {
											sumBestBlastGoAnnotationsSimpleScore += bestBlastGoAnnotationSimpleScore;
											covBestBlastGoAnnotationsSimpleScore++;
										}
									} else {
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
									}
								}
								averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsSimpleScore/covBestBlastGoAnnotationsSimpleScore);
								coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsSimpleScore/getProteins().size());
								if (getSettings().doWriteFscoreDetailsToOutput()) {
									columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-Simple-F-Score-Precision";
									double sumBestBlastGoAnnotationsSimpleScorePrecision = 0.0;
									int covBestBlastGoAnnotationsSimpleScorePrecision = 0;
									for (Protein prot : getProteins()) {
										if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
											BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
											Double bestBlastGoAnnotationSimpleScorePrecision = bestBr.getSimpleGoAnnotationScore().getPrecision();
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationSimpleScorePrecision));
											if (!bestBlastGoAnnotationSimpleScorePrecision.isNaN()) {
												sumBestBlastGoAnnotationsSimpleScorePrecision += bestBlastGoAnnotationSimpleScorePrecision;
												covBestBlastGoAnnotationsSimpleScorePrecision++;
											}
										} else {
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
										}
									}
									averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsSimpleScorePrecision/covBestBlastGoAnnotationsSimpleScorePrecision);
									coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsSimpleScorePrecision/getProteins().size());
									columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-Simple-F-Score-Recall";
									double sumBestBlastGoAnnotationsSimpleScoreRecall = 0.0;
									int covBestBlastGoAnnotationsSimpleScoreRecall = 0;
									for (Protein prot : getProteins()) {
										if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
											BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
											Double bestBlastGoAnnotationSimpleScoreRecall = bestBr.getSimpleGoAnnotationScore().getRecall();
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationSimpleScoreRecall));
											if (!bestBlastGoAnnotationSimpleScoreRecall.isNaN()) {
												sumBestBlastGoAnnotationsSimpleScoreRecall += bestBlastGoAnnotationSimpleScoreRecall;
												covBestBlastGoAnnotationsSimpleScoreRecall++;
											}
										} else {
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
										}
									}
									averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsSimpleScoreRecall/covBestBlastGoAnnotationsSimpleScoreRecall);
									coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsSimpleScoreRecall/getProteins().size());
								}
							}
							if (getSettings().doCalculateAncestryGoF1Scores()) {
								columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-Ancestry-F-Score";
								double sumBestBlastGoAnnotationsAncestryScore = 0.0;
								int covBestBlastGoAnnotationsAncestryScore = 0;
								for (Protein prot : getProteins()) {
									if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
										BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
										Double bestBlastGoAnnotationAncestryScore = bestBr.getAncestryGoAnnotationScore().getScore();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationAncestryScore));
										if (!bestBlastGoAnnotationAncestryScore.isNaN()) {
											sumBestBlastGoAnnotationsAncestryScore += bestBlastGoAnnotationAncestryScore;
											covBestBlastGoAnnotationsAncestryScore++;
										}
									} else {
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
									}
								}
								averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsAncestryScore/covBestBlastGoAnnotationsAncestryScore);
								coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsAncestryScore/getProteins().size());
								if (getSettings().doWriteFscoreDetailsToOutput()) {
									columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-Ancestry-F-Score-Precision";
									double sumBestBlastGoAnnotationsAncestryScorePrecision = 0.0;
									int covBestBlastGoAnnotationsAncestryScorePrecision = 0;
									for (Protein prot : getProteins()) {
										if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
											BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
											Double bestBlastGoAnnotationAncestryScorePrecision = bestBr.getAncestryGoAnnotationScore().getPrecision();
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationAncestryScorePrecision));
											if (!bestBlastGoAnnotationAncestryScorePrecision.isNaN()) {
												sumBestBlastGoAnnotationsAncestryScorePrecision += bestBlastGoAnnotationAncestryScorePrecision;
												covBestBlastGoAnnotationsAncestryScorePrecision++;
											}
										} else {
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
										}
									}
									averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsAncestryScorePrecision/covBestBlastGoAnnotationsAncestryScorePrecision);
									coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsAncestryScorePrecision/getProteins().size());
									columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-Ancestry-F-Score-Recall";
									double sumBestBlastGoAnnotationsAncestryScoreRecall = 0.0;
									int covBestBlastGoAnnotationsAncestryScoreRecall = 0;
									for (Protein prot : getProteins()) {
										if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
											BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
											Double bestBlastGoAnnotationAncestryScoreRecall = bestBr.getAncestryGoAnnotationScore().getRecall();
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationAncestryScoreRecall));
											if (!bestBlastGoAnnotationAncestryScoreRecall.isNaN()) {
												sumBestBlastGoAnnotationsAncestryScoreRecall += bestBlastGoAnnotationAncestryScoreRecall;
												covBestBlastGoAnnotationsAncestryScoreRecall++;
											}
										} else {
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
										}
									}
									averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsAncestryScoreRecall/covBestBlastGoAnnotationsAncestryScoreRecall);
									coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsAncestryScoreRecall/getProteins().size());
								}
							}
							if (getSettings().doCalculateSemSimGoF1Scores()) {
								columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-SemSim-F-Score";
								double sumBestBlastGoAnnotationsSemSimScore = 0.0;
								int covBestBlastGoAnnotationsSemSimScore = 0;
								for (Protein prot : getProteins()) {
									if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
										BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
										Double bestBlastGoAnnotationSemSimScore = bestBr.getSemSimGoAnnotationScore().getScore();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationSemSimScore));
										if (!bestBlastGoAnnotationSemSimScore.isNaN()) {
											sumBestBlastGoAnnotationsSemSimScore += bestBlastGoAnnotationSemSimScore;
											covBestBlastGoAnnotationsSemSimScore++;
										}
									} else {
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
									}
								}
								averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsSemSimScore/covBestBlastGoAnnotationsSemSimScore);
								coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsSemSimScore/getProteins().size());
								if (getSettings().doWriteFscoreDetailsToOutput()) {
									columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-SemSim-F-Score-Precision";
									double sumBestBlastGoAnnotationsSemSimScorePrecision = 0.0;
									int covBestBlastGoAnnotationsSemSimScorePrecision = 0;
									for (Protein prot : getProteins()) {
										if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
											BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
											Double bestBlastGoAnnotationSemSimScorePrecision = bestBr.getSemSimGoAnnotationScore().getPrecision();
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationSemSimScorePrecision));
											if (!bestBlastGoAnnotationSemSimScorePrecision.isNaN()) {
												sumBestBlastGoAnnotationsSemSimScorePrecision += bestBlastGoAnnotationSemSimScorePrecision;
												covBestBlastGoAnnotationsSemSimScorePrecision++;
											}
										} else {
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
										}
									}
									averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsSemSimScorePrecision/covBestBlastGoAnnotationsSemSimScorePrecision);
									coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsSemSimScorePrecision/getProteins().size());
									columnNamesLine += separator + "Best-BlastHit-against-" + blastDb + "-GO-Annotations-SemSim-F-Score-Recall";
									double sumBestBlastGoAnnotationsSemSimScoreRecall = 0.0;
									int covBestBlastGoAnnotationsSemSimScoreRecall = 0;
									for (Protein prot : getProteins()) {
										if (prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb) != null) {
											BlastResult bestBr = prot.getEvaluationScoreCalculator().getBestUnchangedBlastResults().get(blastDb);
											Double bestBlastGoAnnotationSemSimScoreRecall = bestBr.getSemSimGoAnnotationScore().getRecall();
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(bestBlastGoAnnotationSemSimScoreRecall));
											if (!bestBlastGoAnnotationSemSimScoreRecall.isNaN()) {
												sumBestBlastGoAnnotationsSemSimScoreRecall += bestBlastGoAnnotationSemSimScoreRecall;
												covBestBlastGoAnnotationsSemSimScoreRecall++;
											}
										} else {
											proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
										}
									}
									averagesLine += separator + FRMT.format(sumBestBlastGoAnnotationsSemSimScoreRecall/covBestBlastGoAnnotationsSemSimScoreRecall);
									coveragesLine += separator + FRMT.format((double)covBestBlastGoAnnotationsSemSimScoreRecall/getProteins().size());
								}
							}
						}
					}
				}
			}
		}
		
		// Competitors
		if (getSettings().hasCompetitors()) {
			List<String> competitorList = new ArrayList<>(getSettings().getCompetitorSettings().keySet());
			Collections.sort(competitorList);
			for (String competitor : competitorList) {
				// Competitor description evaluation 
				columnNamesLine += separator + competitor + "-Description" + separator + competitor + "-Description-Length" + separator + competitor + "-Description-Evaluation-Score";
				int sumCompetitorDescriptionLength = 0;
				int covCompetitorDescriptionLength = 0;
				double sumCompetitorDescriptionScore = 0.0;
				int covCompetitorDescriptionScore = 0;
				for (Protein prot : getProteins()) {
					Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
					if (compAnnots != null) {
						CompetitorAnnotation annot = compAnnots.get(competitor);
						if (annot != null) {
							int competitorDescriptionLength = annot.getEvaluationTokens().size();
							Double competitorDescriptionScore = annot.getEvaluationScore().getScore();
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + annot.getDescription() + separator + competitorDescriptionLength + separator + FRMT.format(competitorDescriptionScore));
							if (competitorDescriptionLength > 0) {
								sumCompetitorDescriptionLength += annot.getEvaluationTokens().size();
								covCompetitorDescriptionLength++;
							}
							if (!competitorDescriptionScore.isNaN()) {
								sumCompetitorDescriptionScore += competitorDescriptionScore;
								covCompetitorDescriptionScore++;
							}
						} else {
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator + separator);
						}
					} else {
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator + separator);
					}
				}
				averagesLine += separator + separator + FRMT.format((double)sumCompetitorDescriptionLength/covCompetitorDescriptionLength) + separator + FRMT.format(sumCompetitorDescriptionScore/covCompetitorDescriptionScore);
				coveragesLine += separator + separator + FRMT.format((double)covCompetitorDescriptionLength/getProteins().size()) + separator + FRMT.format((double)covCompetitorDescriptionScore/getProteins().size());
				// Competitor description evaluation details
				if (getSettings().doWriteFscoreDetailsToOutput()) {
					columnNamesLine += separator + competitor + "-Description-Evaluation-Score-Precision";
					columnNamesLine += separator + competitor + "-Description-Evaluation-Score-Recall";
					double sumCompetitorDescriptionScorePrecision = 0.0;
					int covCompetitorDescriptionScorePrecision = 0;
					double sumCompetitorDescriptionScoreRecall = 0.0;
					int covCompetitorDescriptionScoreRecall = 0;
					for (Protein prot : getProteins()) {
						Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
						if (compAnnots != null) {
							CompetitorAnnotation annot = compAnnots.get(competitor);
							if (annot != null) {
								Double CompetitorDescriptionScorePrecision = annot.getEvaluationScore().getPrecision();
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(CompetitorDescriptionScorePrecision));
								if (!CompetitorDescriptionScorePrecision.isNaN()) {
									sumCompetitorDescriptionScorePrecision += CompetitorDescriptionScorePrecision;
									covCompetitorDescriptionScorePrecision++;
								}
								Double CompetitorDescriptionScoreRecall = annot.getEvaluationScore().getRecall();
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(CompetitorDescriptionScoreRecall));
								if (!CompetitorDescriptionScoreRecall.isNaN()) {
									sumCompetitorDescriptionScoreRecall += CompetitorDescriptionScoreRecall;
									covCompetitorDescriptionScoreRecall++;
								}
							} else {
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
							}
						} else {
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
						}
					}
					averagesLine += separator + FRMT.format(sumCompetitorDescriptionScorePrecision/covCompetitorDescriptionScorePrecision) + separator + FRMT.format(sumCompetitorDescriptionScoreRecall/covCompetitorDescriptionScoreRecall);
					coveragesLine += separator + FRMT.format((double)covCompetitorDescriptionScorePrecision/getProteins().size()) + separator + FRMT.format((double)covCompetitorDescriptionScoreRecall/getProteins().size());
				}
				// Competitor GO terms
				if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
					columnNamesLine += separator + competitor + "-GO-Annotations";
					for (Protein prot : getProteins()) {
						Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
						if (compAnnots != null) {
							CompetitorAnnotation annot = compAnnots.get(competitor);
							if (annot != null) {
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + combineGoTermsToString(annot.getGoAnnotations()));
							} else {
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
							}
						} else {
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
						}
					}
					averagesLine += separator;
					coveragesLine += separator;
					// Competitor GO term evaluation
					if (getSettings().doCalculateSimpleGoF1Scores()){
						columnNamesLine += separator + competitor + "-GO-Annotations-Simple-F-Score";
						double sumCompetitorSimpleGoScore = 0.0;
						int covCompetitorSimpleGoScore = 0;
						for (Protein prot : getProteins()) {
							Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
							if (compAnnots != null) {
								CompetitorAnnotation annot = compAnnots.get(competitor);
								if (annot != null) {
									Double competitorSimpleGoScore = annot.getSimpleGoAnnotationScore().getScore();
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorSimpleGoScore));
									if (!competitorSimpleGoScore.isNaN()) {
										sumCompetitorSimpleGoScore += competitorSimpleGoScore;
										covCompetitorSimpleGoScore++;
									}
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
								}
							} else {
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
							}
						}
						averagesLine += separator + FRMT.format(sumCompetitorSimpleGoScore/covCompetitorSimpleGoScore);
						coveragesLine += separator + FRMT.format((double)covCompetitorSimpleGoScore/getProteins().size());
						if (getSettings().doWriteFscoreDetailsToOutput()) {
							columnNamesLine += separator + competitor + "-GO-Annotations-Simple-F-Score-Precision";
							columnNamesLine += separator + competitor + "-GO-Annotations-Simple-F-Score-Recall";
							double sumCompetitorSimpleGoScorePrecision = 0.0;
							int covCompetitorSimpleGoScorePrecision = 0;
							double sumCompetitorSimpleGoScoreRecall = 0.0;
							int covCompetitorSimpleGoScoreRecall = 0;
							for (Protein prot : getProteins()) {
								Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
								if (compAnnots != null) {
									CompetitorAnnotation annot = compAnnots.get(competitor);
									if (annot != null) {
										Double competitorSimpleGoScorePrecision = annot.getSimpleGoAnnotationScore().getPrecision();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorSimpleGoScorePrecision));
										if (!competitorSimpleGoScorePrecision.isNaN()) {
											sumCompetitorSimpleGoScorePrecision += competitorSimpleGoScorePrecision;
											covCompetitorSimpleGoScorePrecision++;
										}
										Double competitorSimpleGoScoreRecall = annot.getSimpleGoAnnotationScore().getRecall();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorSimpleGoScoreRecall));
										if (!competitorSimpleGoScoreRecall.isNaN()) {
											sumCompetitorSimpleGoScoreRecall += competitorSimpleGoScoreRecall;
											covCompetitorSimpleGoScoreRecall++;
										}
									} else {
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
									}
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
								}
							}
							averagesLine += separator + FRMT.format(sumCompetitorSimpleGoScorePrecision/covCompetitorSimpleGoScorePrecision) + separator + FRMT.format(sumCompetitorSimpleGoScoreRecall/covCompetitorSimpleGoScoreRecall);
							coveragesLine += separator + FRMT.format((double)covCompetitorSimpleGoScorePrecision/getProteins().size()) + separator + FRMT.format((double)covCompetitorSimpleGoScoreRecall/getProteins().size());
						}
					}
					if (getSettings().doCalculateAncestryGoF1Scores()){
						columnNamesLine += separator + competitor + "-GO-Annotations-Ancestry-F-Score";
						double sumCompetitorAncestryGoScore = 0.0;
						int covCompetitorAncestryGoScore = 0;
						for (Protein prot : getProteins()) {
							Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
							if (compAnnots != null) {
								CompetitorAnnotation annot = compAnnots.get(competitor);
								if (annot != null) {
									Double competitorAncestryGoScore = annot.getAncestryGoAnnotationScore().getScore();
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorAncestryGoScore));
									if (!competitorAncestryGoScore.isNaN()) {
										sumCompetitorAncestryGoScore += competitorAncestryGoScore;
										covCompetitorAncestryGoScore++;
									}
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
								}
							} else {
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
							}
						}
						averagesLine += separator + FRMT.format(sumCompetitorAncestryGoScore/covCompetitorAncestryGoScore);
						coveragesLine += separator + FRMT.format((double)covCompetitorAncestryGoScore/getProteins().size());
						if (getSettings().doWriteFscoreDetailsToOutput()) {
							columnNamesLine += separator + competitor + "-GO-Annotations-Ancestry-F-Score-Precision";
							columnNamesLine += separator + competitor + "-GO-Annotations-Ancestry-F-Score-Recall";
							double sumCompetitorAncestryGoScorePrecision = 0.0;
							int covCompetitorAncestryGoScorePrecision = 0;
							double sumCompetitorAncestryGoScoreRecall = 0.0;
							int covCompetitorAncestryGoScoreRecall = 0;
							for (Protein prot : getProteins()) {
								Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
								if (compAnnots != null) {
									CompetitorAnnotation annot = compAnnots.get(competitor);
									if (annot != null) {
										Double competitorAncestryGoScorePrecision = annot.getAncestryGoAnnotationScore().getPrecision();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorAncestryGoScorePrecision));
										if (!competitorAncestryGoScorePrecision.isNaN()) {
											sumCompetitorAncestryGoScorePrecision += competitorAncestryGoScorePrecision;
											covCompetitorAncestryGoScorePrecision++;
										}
										Double competitorAncestryGoScoreRecall = annot.getAncestryGoAnnotationScore().getRecall();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorAncestryGoScoreRecall));
										if (!competitorAncestryGoScoreRecall.isNaN()) {
											sumCompetitorAncestryGoScoreRecall += competitorAncestryGoScoreRecall;
											covCompetitorAncestryGoScoreRecall++;
										}
									} else {
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
									}
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
								}
							}
							averagesLine += separator + FRMT.format(sumCompetitorAncestryGoScorePrecision/covCompetitorAncestryGoScorePrecision) + separator + FRMT.format(sumCompetitorAncestryGoScoreRecall/covCompetitorAncestryGoScoreRecall);
							coveragesLine += separator + FRMT.format((double)covCompetitorAncestryGoScorePrecision/getProteins().size()) + separator + FRMT.format((double)covCompetitorAncestryGoScoreRecall/getProteins().size());
						}
					}
					if (getSettings().doCalculateSemSimGoF1Scores()){
						columnNamesLine += separator + competitor + "-GO-Annotations-SemSim-F-Score";
						double sumCompetitorSemSimGoScore = 0.0;
						int covCompetitorSemSimGoScore = 0;
						for (Protein prot : getProteins()) {
							Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
							if (compAnnots != null) {
								CompetitorAnnotation annot = compAnnots.get(competitor);
								if (annot != null) {
									Double competitorSemSimGoScore = annot.getSemSimGoAnnotationScore().getScore();
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorSemSimGoScore));
									if (!competitorSemSimGoScore.isNaN()) {
										sumCompetitorSemSimGoScore += competitorSemSimGoScore;
										covCompetitorSemSimGoScore++;
									}
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
								}
							} else {
								proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator);
							}
						}
						averagesLine += separator + FRMT.format(sumCompetitorSemSimGoScore/covCompetitorSemSimGoScore);
						coveragesLine += separator + FRMT.format((double)covCompetitorSemSimGoScore/getProteins().size());
						if (getSettings().doWriteFscoreDetailsToOutput()) {
							columnNamesLine += separator + competitor + "-GO-Annotations-SemSim-F-Score-Precision";
							columnNamesLine += separator + competitor + "-GO-Annotations-SemSim-F-Score-Recall";
							double sumCompetitorSemSimGoScorePrecision = 0.0;
							int covCompetitorSemSimGoScorePrecision = 0;
							double sumCompetitorSemSimGoScoreRecall = 0.0;
							int covCompetitorSemSimGoScoreRecall = 0;
							for (Protein prot : getProteins()) {
								Map<String, CompetitorAnnotation> compAnnots = prot.getEvaluationScoreCalculator().getCompetitorAnnotations();
								if (compAnnots != null) {
									CompetitorAnnotation annot = compAnnots.get(competitor);
									if (annot != null) {
										Double competitorSemSimGoScorePrecision = annot.getSemSimGoAnnotationScore().getPrecision();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorSemSimGoScorePrecision));
										if (!competitorSemSimGoScorePrecision.isNaN()) {
											sumCompetitorSemSimGoScorePrecision += competitorSemSimGoScorePrecision;
											covCompetitorSemSimGoScorePrecision++;
										}
										Double competitorSemSimGoScoreRecall = annot.getSemSimGoAnnotationScore().getRecall();
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(competitorSemSimGoScoreRecall));
										if (!competitorSemSimGoScoreRecall.isNaN()) {
											sumCompetitorSemSimGoScoreRecall += competitorSemSimGoScoreRecall;
											covCompetitorSemSimGoScoreRecall++;
										}
									} else {
										proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
									}
								} else {
									proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator);
								}
							}
							averagesLine += separator + FRMT.format(sumCompetitorSemSimGoScorePrecision/covCompetitorSemSimGoScorePrecision) + separator + FRMT.format(sumCompetitorSemSimGoScoreRecall/covCompetitorSemSimGoScoreRecall);
							coveragesLine += separator + FRMT.format((double)covCompetitorSemSimGoScorePrecision/getProteins().size()) + separator + FRMT.format((double)covCompetitorSemSimGoScoreRecall/getProteins().size());
						}
					}
				}
			}
		}
		// Highest possible description evaluation scores
		if (getSettings().doFindHighestPossibleEvaluationScore()) {
			columnNamesLine += separator + "Blast-Result-With-Highest-Possible-Description-Score" + separator + "Highest-Possible-Description-Score";
			double sumHighestPossibleDescriptionScore = 0.0;
			int covHighestPossibleDescriptionScore = 0;
			for (Protein prot : getProteins()) {
				BlastResult br = prot.getEvaluationScoreCalculator().getBlastResultWithHighestPossibleDescriptionScore();
				if (br != null) {
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + "\"" + br.getAccession() + " " + br.getDescription() + "\"");
					Double highestPossibleDescriptionScore = prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionScore().getScore();
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleDescriptionScore));
					if (!highestPossibleDescriptionScore.isNaN()) {
						sumHighestPossibleDescriptionScore += highestPossibleDescriptionScore;
						covHighestPossibleDescriptionScore++;
					}
				} else {
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + separator + "NaN");
				}
			}
			averagesLine += separator + separator + FRMT.format(sumHighestPossibleDescriptionScore/covHighestPossibleDescriptionScore);
			coveragesLine += separator + separator + FRMT.format((double)covHighestPossibleDescriptionScore/getProteins().size());
			if (getSettings().doWriteFscoreDetailsToOutput()) {
				columnNamesLine += separator + "Highest-Possible-Description-Score-Precision" + separator + "Highest-Possible-Description-Score-Recall";
				double sumHighestPossibleDescriptionScorePrecision = 0.0;
				int covHighestPossibleDescriptionScorePrecision = 0;
				double sumHighestPossibleDescriptionScoreRecall = 0.0;
				int covHighestPossibleDescriptionScoreRecall = 0;
				for (Protein prot : getProteins()) {
					BlastResult br = prot.getEvaluationScoreCalculator().getBlastResultWithHighestPossibleDescriptionScore();
					if (br != null) {
						Double highestPossibleDescriptionScorePrecision = prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionScore().getPrecision();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleDescriptionScorePrecision));
						if (!highestPossibleDescriptionScorePrecision.isNaN()) {
							sumHighestPossibleDescriptionScorePrecision+=highestPossibleDescriptionScorePrecision;
							covHighestPossibleDescriptionScorePrecision++;
						}
						Double highestPossibleDescriptionScoreRecall = prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionScore().getRecall();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleDescriptionScoreRecall));
						if (!highestPossibleDescriptionScoreRecall.isNaN()) {
							sumHighestPossibleDescriptionScoreRecall+=highestPossibleDescriptionScoreRecall;
							covHighestPossibleDescriptionScoreRecall++;
						}
					} else {
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + "NaN" + separator + "NaN");
					}
				}
				averagesLine += separator + FRMT.format(sumHighestPossibleDescriptionScorePrecision/covHighestPossibleDescriptionScorePrecision) + separator + FRMT.format(sumHighestPossibleDescriptionScoreRecall/covHighestPossibleDescriptionScoreRecall);
				coveragesLine += separator + FRMT.format((double)covHighestPossibleDescriptionScorePrecision/getProteins().size()) + separator + FRMT.format((double)covHighestPossibleDescriptionScoreRecall/getProteins().size()); 
			}
		}
		// Highest possible description precision
		if (getSettings().doFindHighestPossiblePrecision()) {
			columnNamesLine += separator + "Highest-Possible-Description-Precision";
			double sumHighestPossibleDescriptionPrecision = 0.0;
			int covHighestPossibleDescriptionPrecision = 0;
			for(Protein prot : getProteins()) {
				Double highestPossibleDescriptionPrecision = prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionPrecision();
				proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleDescriptionPrecision));
				if (!highestPossibleDescriptionPrecision.isNaN()) {
					sumHighestPossibleDescriptionPrecision += highestPossibleDescriptionPrecision;
					covHighestPossibleDescriptionPrecision++;
				}
			}
			averagesLine += separator + FRMT.format(sumHighestPossibleDescriptionPrecision/covHighestPossibleDescriptionPrecision);
			coveragesLine += separator + FRMT.format((double)covHighestPossibleDescriptionPrecision/getProteins().size());
		}
		// Highest possible description recall
		if (getSettings().doFindHighestPossibleRecall()) {
			columnNamesLine += separator + "Highest-Possible-Description-Recall";
			double sumHighestPossibleDescriptionRecall = 0.0;
			int covHighestPossibleDescriptionRecall = 0;
			for(Protein prot : getProteins()) {
				Double highestPossibleDescriptionRecall = prot.getEvaluationScoreCalculator().getHighestPossibleDescriptionRecall();
				proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleDescriptionRecall));
				if (!highestPossibleDescriptionRecall.isNaN()) {
					sumHighestPossibleDescriptionRecall += highestPossibleDescriptionRecall;
					covHighestPossibleDescriptionRecall++;
				}
			}
			averagesLine += separator + FRMT.format(sumHighestPossibleDescriptionRecall/covHighestPossibleDescriptionRecall);
			coveragesLine += separator + FRMT.format((double)covHighestPossibleDescriptionRecall/getProteins().size());
		}
		// AHRD's GO term evaluation
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			columnNamesLine += separator + "Ground-Truth-GO-Annotations";
			for (Protein prot : getProteins()) {
				proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + combineGoTermsToString(prot.getEvaluationScoreCalculator().getGroundTruthGoAnnoatations()));
			}
			averagesLine += separator;
			coveragesLine += separator;
			if (getSettings().doCalculateSimpleGoF1Scores()) {
				columnNamesLine += separator + "AHRD-GO-Annotations-Simple-F-Score";
				double sumSimpleGoAnnotationScore = 0.0;
				int covSimpleGoAnnotationScore = 0;
				for (Protein prot : getProteins()) {
					Double simpleGoAnnotationScore = prot.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getScore();
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(simpleGoAnnotationScore));
					if (!simpleGoAnnotationScore.isNaN()) {
						sumSimpleGoAnnotationScore += simpleGoAnnotationScore;
						covSimpleGoAnnotationScore++;
					}
				}
				averagesLine += separator + FRMT.format(sumSimpleGoAnnotationScore/covSimpleGoAnnotationScore);
				coveragesLine += separator + FRMT.format((double)covSimpleGoAnnotationScore/getProteins().size());
				if (getSettings().doWriteFscoreDetailsToOutput()) {
					columnNamesLine += separator + "AHRD-GO-Annotations-Simple-F-Score-Precision" + separator + "AHRD-GO-Annotations-Simple-F-Score-Recall";
					double sumSimpleGoAnnotationScorePrecision = 0.0;
					int covSimpleGoAnnotationScorePrecision = 0;
					double sumSimpleGoAnnotationScoreRecall = 0.0;
					int covSimpleGoAnnotationScoreRecall = 0;
					for (Protein prot : getProteins()) {
						Double simpleGoAnnotationScorePrecision = prot.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getPrecision();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(simpleGoAnnotationScorePrecision));
						if (!simpleGoAnnotationScorePrecision.isNaN()) {
							sumSimpleGoAnnotationScorePrecision += simpleGoAnnotationScorePrecision;
							covSimpleGoAnnotationScorePrecision++;
						}
						Double simpleGoAnnotationScoreRecall = prot.getEvaluationScoreCalculator().getSimpleGoAnnotationScore().getRecall();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(simpleGoAnnotationScoreRecall));
						if (!simpleGoAnnotationScoreRecall.isNaN()) {
							sumSimpleGoAnnotationScoreRecall += simpleGoAnnotationScoreRecall;
							covSimpleGoAnnotationScoreRecall++;
						}
					}
					averagesLine += separator + FRMT.format(sumSimpleGoAnnotationScorePrecision/covSimpleGoAnnotationScorePrecision) + separator + FRMT.format(sumSimpleGoAnnotationScoreRecall/covSimpleGoAnnotationScoreRecall);
					coveragesLine += separator + FRMT.format((double)covSimpleGoAnnotationScorePrecision/getProteins().size()) + separator + FRMT.format((double)covSimpleGoAnnotationScoreRecall/getProteins().size());
				}
			}
			if (getSettings().doCalculateAncestryGoF1Scores()) {
				columnNamesLine += separator + "AHRD-GO-Annotations-Ancestry-F-Score";
				double sumAncestryGoAnnotationScore = 0.0;
				int covAncestryGoAnnotationScore = 0;
				for (Protein prot : getProteins()) {
					Double ancestryGoAnnotationScore = prot.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getScore();
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(ancestryGoAnnotationScore));
					if (!ancestryGoAnnotationScore.isNaN()) {
						sumAncestryGoAnnotationScore += ancestryGoAnnotationScore;
						covAncestryGoAnnotationScore++;
					}
				}
				averagesLine += separator + FRMT.format(sumAncestryGoAnnotationScore/covAncestryGoAnnotationScore);
				coveragesLine += separator + FRMT.format((double)covAncestryGoAnnotationScore/getProteins().size());
				if (getSettings().doWriteFscoreDetailsToOutput()) {
					columnNamesLine += separator + "AHRD-GO-Annotations-Ancestry-F-Score-Precision" + separator + "AHRD-GO-Annotations-Ancestry-F-Score-Recall";
					double sumAncestryGoAnnotationScorePrecision = 0.0;
					int covAncestryGoAnnotationScorePrecision = 0;
					double sumAncestryGoAnnotationScoreRecall = 0.0;
					int covAncestryGoAnnotationScoreRecall = 0;
					for (Protein prot : getProteins()) {
						Double ancestryGoAnnotationScorePrecision = prot.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getPrecision();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(ancestryGoAnnotationScorePrecision));
						if (!ancestryGoAnnotationScorePrecision.isNaN()) {
							sumAncestryGoAnnotationScorePrecision += ancestryGoAnnotationScorePrecision;
							covAncestryGoAnnotationScorePrecision++;
						}
						Double ancestryGoAnnotationScoreRecall = prot.getEvaluationScoreCalculator().getAncestryGoAnnotationScore().getRecall();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(ancestryGoAnnotationScoreRecall));
						if (!ancestryGoAnnotationScoreRecall.isNaN()) {
							sumAncestryGoAnnotationScoreRecall += ancestryGoAnnotationScoreRecall;
							covAncestryGoAnnotationScoreRecall++;
						}
					}
					averagesLine += separator + FRMT.format(sumAncestryGoAnnotationScorePrecision/covAncestryGoAnnotationScorePrecision) + separator + FRMT.format(sumAncestryGoAnnotationScoreRecall/covAncestryGoAnnotationScoreRecall);
					coveragesLine += separator + FRMT.format((double)covAncestryGoAnnotationScorePrecision/getProteins().size()) + separator + FRMT.format((double)covAncestryGoAnnotationScoreRecall/getProteins().size());
				}
			}
			if (getSettings().doCalculateSemSimGoF1Scores()) {
				columnNamesLine += separator + "AHRD-GO-Annotations-SemSim-F-Score";
				double sumSemSimGoAnnotationScore = 0.0;
				int covSemSimGoAnnotationScore = 0;
				for (Protein prot : getProteins()) {
					Double semSimGoAnnotationScore = prot.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getScore();
					proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(semSimGoAnnotationScore));
					if (!semSimGoAnnotationScore.isNaN()) {
						sumSemSimGoAnnotationScore += semSimGoAnnotationScore;
						covSemSimGoAnnotationScore++;
					}
				}
				averagesLine += separator + FRMT.format(sumSemSimGoAnnotationScore/covSemSimGoAnnotationScore);
				coveragesLine += separator + FRMT.format((double)covSemSimGoAnnotationScore/getProteins().size());
				if (getSettings().doWriteFscoreDetailsToOutput()) {
					columnNamesLine += separator + "AHRD-GO-Annotations-SemSim-F-Score-Precision" + separator + "AHRD-GO-Annotations-SemSim-F-Score-Recall";
					double sumSemSimGoAnnotationScorePrecision = 0.0;
					int covSemSimGoAnnotationScorePrecision = 0;
					double sumSemSimGoAnnotationScoreRecall = 0.0;
					int covSemSimGoAnnotationScoreRecall = 0;
					for (Protein prot : getProteins()) {
						Double semSimGoAnnotationScorePrecision = prot.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getPrecision();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(semSimGoAnnotationScorePrecision));
						if (!semSimGoAnnotationScorePrecision.isNaN()) {
							sumSemSimGoAnnotationScorePrecision += semSimGoAnnotationScorePrecision;
							covSemSimGoAnnotationScorePrecision++;
						}
						Double semSimGoAnnotationScoreRecall = prot.getEvaluationScoreCalculator().getSemSimGoAnnotationScore().getRecall();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(semSimGoAnnotationScoreRecall));
						if (!semSimGoAnnotationScoreRecall.isNaN()) {
							sumSemSimGoAnnotationScoreRecall += semSimGoAnnotationScoreRecall;
							covSemSimGoAnnotationScoreRecall++;
						}
					}
					averagesLine += separator + FRMT.format(sumSemSimGoAnnotationScorePrecision/covSemSimGoAnnotationScorePrecision) + separator + FRMT.format(sumSemSimGoAnnotationScoreRecall/covSemSimGoAnnotationScoreRecall);
					coveragesLine += separator + FRMT.format((double)covSemSimGoAnnotationScorePrecision/getProteins().size()) + separator + FRMT.format((double)covSemSimGoAnnotationScoreRecall/getProteins().size());
				}
			}
			// Highest possible go annotation scores
			if (getSettings().doFindHighestPossibleGoScore()) {
				if (getSettings().doCalculateSimpleGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-BlastResult-GO-Annotations-Simple-F-Score";
					double sumHighestPossibleSimpleGoAnnotationScore = 0.0;
					int covHighestPossibleSimpleGoAnnotationScore = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleSimpleGoAnnotationScore = prot.getEvaluationScoreCalculator().getHighestPossibleSimpleGoAnnotationScore().getScore();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSimpleGoAnnotationScore));
						if (!highestPossibleSimpleGoAnnotationScore.isNaN()) {
							sumHighestPossibleSimpleGoAnnotationScore += highestPossibleSimpleGoAnnotationScore;
							covHighestPossibleSimpleGoAnnotationScore++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleSimpleGoAnnotationScore/covHighestPossibleSimpleGoAnnotationScore);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleSimpleGoAnnotationScore/getProteins().size());
					if (getSettings().doWriteFscoreDetailsToOutput()) {
						columnNamesLine += separator + "Highest-Possible-BlastResult-GO-Annotations-Simple-F-Score-Precision" + separator + "Highest-Possible-BlastResult-GO-Annotations-Simple-F-Score-Recall";
						double sumHighestPossibleSimpleGoAnnotationScorePrecision = 0.0;
						int covHighestPossibleSimpleGoAnnotationScorePrecision = 0;
						double sumHighestPossibleSimpleGoAnnotationScoreRecall = 0.0;
						int covHighestPossibleSimpleGoAnnotationScoreRecall = 0;
						for (Protein prot : getProteins()) {
							Double highestPossibleSimpleGoAnnotationScorePrecision = prot.getEvaluationScoreCalculator().getHighestPossibleSimpleGoAnnotationScore().getPrecision();
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSimpleGoAnnotationScorePrecision));
							if (!highestPossibleSimpleGoAnnotationScorePrecision.isNaN()) {
								sumHighestPossibleSimpleGoAnnotationScorePrecision += highestPossibleSimpleGoAnnotationScorePrecision;
								covHighestPossibleSimpleGoAnnotationScorePrecision++;
							}
							Double highestPossibleSimpleGoAnnotationScoreRecall = prot.getEvaluationScoreCalculator().getHighestPossibleSimpleGoAnnotationScore().getRecall();
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSimpleGoAnnotationScoreRecall));
							if (!highestPossibleSimpleGoAnnotationScoreRecall.isNaN()) {
								sumHighestPossibleSimpleGoAnnotationScoreRecall += highestPossibleSimpleGoAnnotationScoreRecall;
								covHighestPossibleSimpleGoAnnotationScoreRecall++;
							}
						}
						averagesLine += separator + FRMT.format(sumHighestPossibleSimpleGoAnnotationScorePrecision/covHighestPossibleSimpleGoAnnotationScorePrecision) + separator + FRMT.format(sumHighestPossibleSimpleGoAnnotationScoreRecall/covHighestPossibleSimpleGoAnnotationScoreRecall);
						coveragesLine += separator + FRMT.format((double)covHighestPossibleSimpleGoAnnotationScorePrecision/getProteins().size()) + separator + FRMT.format((double)covHighestPossibleSimpleGoAnnotationScoreRecall/getProteins().size());
					}
				}
				if (getSettings().doCalculateAncestryGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-BlastResult-GO-Annotations-Ancestry-F-Score";
					double sumHighestPossibleAncestryGoAnnotationScore = 0.0;
					int covHighestPossibleAncestryGoAnnotationScore = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleAncestryGoAnnotationScore = prot.getEvaluationScoreCalculator().getHighestPossibleAncestryGoAnnotationScore().getScore();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleAncestryGoAnnotationScore));
						if (!highestPossibleAncestryGoAnnotationScore.isNaN()) {
							sumHighestPossibleAncestryGoAnnotationScore += highestPossibleAncestryGoAnnotationScore;
							covHighestPossibleAncestryGoAnnotationScore++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleAncestryGoAnnotationScore/covHighestPossibleAncestryGoAnnotationScore);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleAncestryGoAnnotationScore/getProteins().size());
					if (getSettings().doWriteFscoreDetailsToOutput()) {
						columnNamesLine += separator + "Highest-Possible-BlastResult-GO-Annotations-Ancestry-F-Score-Precision" + separator + "Highest-Possible-BlastResult-GO-Annotations-Ancestry-F-Score-Recall";
						double sumHighestPossibleAncestryGoAnnotationScorePrecision = 0.0;
						int covHighestPossibleAncestryGoAnnotationScorePrecision = 0;
						double sumHighestPossibleAncestryGoAnnotationScoreRecall = 0.0;
						int covHighestPossibleAncestryGoAnnotationScoreRecall = 0;
						for (Protein prot : getProteins()) {
							Double highestPossibleAncestryGoAnnotationScorePrecision = prot.getEvaluationScoreCalculator().getHighestPossibleAncestryGoAnnotationScore().getPrecision();
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleAncestryGoAnnotationScorePrecision));
							if (!highestPossibleAncestryGoAnnotationScorePrecision.isNaN()) {
								sumHighestPossibleAncestryGoAnnotationScorePrecision += highestPossibleAncestryGoAnnotationScorePrecision;
								covHighestPossibleAncestryGoAnnotationScorePrecision++;
							}
							Double highestPossibleAncestryGoAnnotationScoreRecall = prot.getEvaluationScoreCalculator().getHighestPossibleAncestryGoAnnotationScore().getRecall();
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleAncestryGoAnnotationScoreRecall));
							if (!highestPossibleAncestryGoAnnotationScoreRecall.isNaN()) {
								sumHighestPossibleAncestryGoAnnotationScoreRecall += highestPossibleAncestryGoAnnotationScoreRecall;
								covHighestPossibleAncestryGoAnnotationScoreRecall++;
							}
						}
						averagesLine += separator + FRMT.format(sumHighestPossibleAncestryGoAnnotationScorePrecision/covHighestPossibleAncestryGoAnnotationScorePrecision) + separator + FRMT.format(sumHighestPossibleAncestryGoAnnotationScoreRecall/covHighestPossibleAncestryGoAnnotationScoreRecall);
						coveragesLine += separator + FRMT.format((double)covHighestPossibleAncestryGoAnnotationScorePrecision/getProteins().size()) + separator + FRMT.format((double)covHighestPossibleAncestryGoAnnotationScoreRecall/getProteins().size());
					}
				}
				if (getSettings().doCalculateSemSimGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-BlastResult-GO-Annotations-SemSim-F-Score";
					double sumHighestPossibleSemSimGoAnnotationScore = 0.0;
					int covHighestPossibleSemSimGoAnnotationScore = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleSemSimGoAnnotationScore = prot.getEvaluationScoreCalculator().getHighestPossibleSemSimGoAnnotationScore().getScore();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSemSimGoAnnotationScore));
						if (!highestPossibleSemSimGoAnnotationScore.isNaN()) {
							sumHighestPossibleSemSimGoAnnotationScore += highestPossibleSemSimGoAnnotationScore;
							covHighestPossibleSemSimGoAnnotationScore++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleSemSimGoAnnotationScore/covHighestPossibleSemSimGoAnnotationScore);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleSemSimGoAnnotationScore/getProteins().size());
					if (getSettings().doWriteFscoreDetailsToOutput()) {
						columnNamesLine += separator + "Highest-Possible-BlastResult-GO-Annotations-SemSim-F-Score-Precision" + separator + "Highest-Possible-BlastResult-GO-Annotations-SemSim-F-Score-Recall";
						double sumHighestPossibleSemSimGoAnnotationScorePrecision = 0.0;
						int covHighestPossibleSemSimGoAnnotationScorePrecision = 0;
						double sumHighestPossibleSemSimGoAnnotationScoreRecall = 0.0;
						int covHighestPossibleSemSimGoAnnotationScoreRecall = 0;
						for (Protein prot : getProteins()) {
							Double highestPossibleSemSimGoAnnotationScorePrecision = prot.getEvaluationScoreCalculator().getHighestPossibleSemSimGoAnnotationScore().getPrecision();
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSemSimGoAnnotationScorePrecision));
							if (!highestPossibleSemSimGoAnnotationScorePrecision.isNaN()) {
								sumHighestPossibleSemSimGoAnnotationScorePrecision += highestPossibleSemSimGoAnnotationScorePrecision;
								covHighestPossibleSemSimGoAnnotationScorePrecision++;
							}
							Double highestPossibleSemSimGoAnnotationScoreRecall = prot.getEvaluationScoreCalculator().getHighestPossibleSemSimGoAnnotationScore().getRecall();
							proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSemSimGoAnnotationScoreRecall));
							if (!highestPossibleSemSimGoAnnotationScoreRecall.isNaN()) {
								sumHighestPossibleSemSimGoAnnotationScoreRecall += highestPossibleSemSimGoAnnotationScoreRecall;
								covHighestPossibleSemSimGoAnnotationScoreRecall++;
							}
						}
						averagesLine += separator + FRMT.format(sumHighestPossibleSemSimGoAnnotationScorePrecision/covHighestPossibleSemSimGoAnnotationScorePrecision) + separator + FRMT.format(sumHighestPossibleSemSimGoAnnotationScoreRecall/covHighestPossibleSemSimGoAnnotationScoreRecall);
						coveragesLine += separator + FRMT.format((double)covHighestPossibleSemSimGoAnnotationScorePrecision/getProteins().size()) + separator + FRMT.format((double)covHighestPossibleSemSimGoAnnotationScoreRecall/getProteins().size());
					}
				}
			}

			// Highest possible Go annotation precision
			if (getSettings().doFindHighestPossiblePrecision()) {
				if (getSettings().doCalculateSimpleGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-Simple-GO-Annotations-Precision";
					double sumHighestPossibleSimpleGoAnnotationPrecision = 0.0;
					int covHighestPossibleSimpleGoAnnotationPrecision = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleSimpleGoAnnotationPrecision = prot.getEvaluationScoreCalculator().getHighestPossibleSimpleGoAnnotationPrecision();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSimpleGoAnnotationPrecision));
						if (!highestPossibleSimpleGoAnnotationPrecision.isNaN()) {
							sumHighestPossibleSimpleGoAnnotationPrecision += highestPossibleSimpleGoAnnotationPrecision;
							covHighestPossibleSimpleGoAnnotationPrecision++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleSimpleGoAnnotationPrecision/covHighestPossibleSimpleGoAnnotationPrecision);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleSimpleGoAnnotationPrecision/getProteins().size());
				}
				if (getSettings().doCalculateAncestryGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-Ancestry-GO-Annotations-Precision";
					double sumHighestPossibleAncestryGoAnnotationPrecision = 0.0;
					int covHighestPossibleAncestryGoAnnotationPrecision = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleAncestryGoAnnotationPrecision = prot.getEvaluationScoreCalculator().getHighestPossibleAncestryGoAnnotationPrecision();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleAncestryGoAnnotationPrecision));
						if (!highestPossibleAncestryGoAnnotationPrecision.isNaN()) {
							sumHighestPossibleAncestryGoAnnotationPrecision += highestPossibleAncestryGoAnnotationPrecision;
							covHighestPossibleAncestryGoAnnotationPrecision++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleAncestryGoAnnotationPrecision/covHighestPossibleAncestryGoAnnotationPrecision);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleAncestryGoAnnotationPrecision/getProteins().size());
				}
				if (getSettings().doCalculateSemSimGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-SemSim-GO-Annotations-Precision";
					double sumHighestPossibleSemSimGoAnnotationPrecision = 0.0;
					int covHighestPossibleSemSimGoAnnotationPrecision = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleSemSimGoAnnotationPrecision = prot.getEvaluationScoreCalculator().getHighestPossibleSemSimGoAnnotationPrecision();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSemSimGoAnnotationPrecision));
						if (!highestPossibleSemSimGoAnnotationPrecision.isNaN()) {
							sumHighestPossibleSemSimGoAnnotationPrecision += highestPossibleSemSimGoAnnotationPrecision;
							covHighestPossibleSemSimGoAnnotationPrecision++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleSemSimGoAnnotationPrecision/covHighestPossibleSemSimGoAnnotationPrecision);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleSemSimGoAnnotationPrecision/getProteins().size());
				}
			}
			if (getSettings().doFindHighestPossibleRecall()) {
				if (getSettings().doCalculateSimpleGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-Simple-GO-Annotations-Recall";
					double sumHighestPossibleSimpleGoAnnotationRecall = 0.0;
					int covHighestPossibleSimpleGoAnnotationRecall = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleSimpleGoAnnotationRecall = prot.getEvaluationScoreCalculator().getHighestPossibleSimpleGoAnnotationRecall();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSimpleGoAnnotationRecall));
						if (!highestPossibleSimpleGoAnnotationRecall.isNaN()) {
							sumHighestPossibleSimpleGoAnnotationRecall += highestPossibleSimpleGoAnnotationRecall;
							covHighestPossibleSimpleGoAnnotationRecall++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleSimpleGoAnnotationRecall/covHighestPossibleSimpleGoAnnotationRecall);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleSimpleGoAnnotationRecall/getProteins().size());
				}
				if (getSettings().doCalculateAncestryGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-Ancestry-GO-Annotations-Recall";
					double sumHighestPossibleAncestryGoAnnotationRecall = 0.0;
					int covHighestPossibleAncestryGoAnnotationRecall = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleAncestryGoAnnotationRecall = prot.getEvaluationScoreCalculator().getHighestPossibleAncestryGoAnnotationRecall();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleAncestryGoAnnotationRecall));
						if (!highestPossibleAncestryGoAnnotationRecall.isNaN()) {
							sumHighestPossibleAncestryGoAnnotationRecall += highestPossibleAncestryGoAnnotationRecall;
							covHighestPossibleAncestryGoAnnotationRecall++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleAncestryGoAnnotationRecall/covHighestPossibleAncestryGoAnnotationRecall);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleAncestryGoAnnotationRecall/getProteins().size());
				}
				if (getSettings().doCalculateSemSimGoF1Scores()) {
					columnNamesLine += separator + "Highest-Possible-SemSim-GO-Annotations-Recall";
					double sumHighestPossibleSemSimGoAnnotationRecall = 0.0;
					int covHighestPossibleSemSimGoAnnotationRecall = 0;
					for (Protein prot : getProteins()) {
						Double highestPossibleSemSimGoAnnotationRecall = prot.getEvaluationScoreCalculator().getHighestPossibleSemSimGoAnnotationRecall();
						proteinOutputLines.replace(prot, proteinOutputLines.get(prot) + separator + FRMT.format(highestPossibleSemSimGoAnnotationRecall));
						if (!highestPossibleSemSimGoAnnotationRecall.isNaN()) {
							sumHighestPossibleSemSimGoAnnotationRecall += highestPossibleSemSimGoAnnotationRecall;
							covHighestPossibleSemSimGoAnnotationRecall++;
						}
					}
					averagesLine += separator + FRMT.format(sumHighestPossibleSemSimGoAnnotationRecall/covHighestPossibleSemSimGoAnnotationRecall);
					coveragesLine += separator + FRMT.format((double)covHighestPossibleSemSimGoAnnotationRecall/getProteins().size());
				}
			}
		}
	}
		
	/**
	 * Creates AHRD's most basic output data
	 */
	private String buildDescription(Protein protein) {
		String descriptionLine = protein.getAccession() + separator;
		// Blast-Results
		if (protein.getDescriptionScoreCalculator().getHighestScoringBlastResult() != null) {
			BlastResult br = protein.getDescriptionScoreCalculator().getHighestScoringBlastResult();
			descriptionLine += br.getAccession() + separator + qualityCode(protein) + separator + br.getDescription();
		} else {
			// Maintain Table's Column-Structure if AHRD can't annotate the proteins 
			descriptionLine += separator + separator + "Unknown protein";
		}
		return descriptionLine;
	}
	
}
