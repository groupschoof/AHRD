package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import ahrd.controller.AHRD;
import ahrd.model.Blast2GoAnnot;
import ahrd.model.BlastResult;
import ahrd.model.GOterm;
import ahrd.model.Protein;

public class EvaluatorOutputWriter extends TsvOutputWriter {

	protected BufferedWriter hrdScoresWriter;

	public EvaluatorOutputWriter(Collection<Protein> proteins) {
		super(proteins);
	}

	public void writeOutput() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(getSettings().getPathToOutput()));
		if (getSettings().doWriteHRDScoresToOutput())
			writeHRDScoresOutputHeader();

		// Header
		bw.write("# AHRD-Version " + AHRD.VERSION + " (Evaluator)\n");
		bw.write("\n");
		// Column-Names:
		bw.write(ahrdColumnNames());
		bw.write("\tHRD-Length\tReference-Description\tRef-Lenght\tEvaluation-Score\tDiff-to-bestCompetitor\tTPR\tFPR");
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			bw.write(buildBestBlastHitsHeader());
		}
		if (getSettings().getWriteTokenSetToOutput()) {
			bw.write("\t\"Tokens (tkn->score)\"");
		}
		if (getSettings().getWriteScoresToOutput()) {
			bw.write("\tSum(Token-Scores)\tTokenHighScore\tCorrection-Factor\tLexical-Score\tRelativeBitScore");
		}
		if (getSettings().hasBlast2GoAnnotations()) {
			bw.write("\tBlast2GO-Description\tBlast2GO-Length\tBlast2GO-Evaluation-Score");
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
				bw.write("\tBlast2GO-Annotations");
				if (getSettings().getCalculateSimpleGoF1Scores())
					bw.write("\tBlast2GO-Annotations-Simple-F-Score");
				if (getSettings().getCalculateAncestryGoF1Scores())
					bw.write("\tBlast2GO-Annotations-Ancestry-F-Score");
				if (getSettings().getCalculateSemSimGoF1Scores())
					bw.write("\tBlast2GO-Annotations-SemSim-F-Score");
			}
		}
		if (getSettings().doFindHighestPossibleEvaluationScore()) {
			bw.write("\tHighest-Blast-Hit-Evaluation-Score");
		}
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
			bw.write("\tReference-GO-Annotations");
			if (getSettings().getCalculateSimpleGoF1Scores())
				bw.write("\tAHRD-GO-Annotations-Simple-F-Score");
			if (getSettings().getCalculateAncestryGoF1Scores())
				bw.write("\tAHRD-GO-Annotations-Ancestry-F-Score");
			if (getSettings().getCalculateSemSimGoF1Scores())
				bw.write("\tAHRD-GO-Annotations-SemSim-F-Score");
		}
		bw.write("\n");

		for (Protein prot : getProteins()) {
			// Generate the Human Readable Description:
			String csvRow = buildDescriptionLine(prot, "\t");

			// Write out the Evaluator-Score and the Reference-Description:
			csvRow += buildTrainerColumns(prot);
			// Append further information, if requested:
			if (getSettings().getWriteBestBlastHitsToOutput()) {
				csvRow += buildBestBlastHitsColumns(prot);
			}
			if (getSettings().getWriteTokenSetToOutput()) {
				csvRow += buildTokenSetCell(prot);
			}
			if (getSettings().getWriteScoresToOutput()) {
				csvRow += buildDescScoreCells(prot);
			}
			if (getSettings().hasBlast2GoAnnotations()) {
				csvRow += buildBlast2GoColumns(prot);
			}
			if (getSettings().doFindHighestPossibleEvaluationScore()) {
				csvRow += buildHighestPossibleEvaluationScoreColumn(prot);
			}
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
				csvRow += buildReferenceGoAnnotationColumns(prot);
			}

			// Write row to CSV:
			csvRow += "\n";
			bw.write(csvRow);

			// If AHRD is requested to write out the AHRD-Score of each
			// BlastHit's Description, do so into another file:
			if (getSettings().doWriteHRDScoresToOutput())
				writeHrdScoresOutput(prot);
		}

		// CLEAN UP:
		bw.close();
		if (getSettings().doWriteHRDScoresToOutput())
			this.hrdScoresWriter.close();
	}

	public String buildHighestPossibleEvaluationScoreColumn(Protein prot) {
		return "\t" + FRMT.format(prot.getEvaluationScoreCalculator().getHighestPossibleEvaluationScore());
	}

	public String buildBlast2GoColumns(Protein prot) {
		String csvCols = "";
		List<Blast2GoAnnot> rankedBlast2GoAnnots = prot.getEvaluationScoreCalculator().sortBlast2GoAnnotsByEvalScore();
		if (rankedBlast2GoAnnots != null && !rankedBlast2GoAnnots.isEmpty()) {
			Blast2GoAnnot bestB2ga = rankedBlast2GoAnnots.get(rankedBlast2GoAnnots.size() - 1);
			csvCols += "\t" + bestB2ga.getDescription() + "\t" + bestB2ga.getEvaluationTokens().size() + "\t"
					+ FRMT.format(bestB2ga.getEvaluationScore());
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
				csvCols += "\t" + combineGoTermsToString(bestB2ga.getGoAnnotations());
				if (getSettings().getCalculateSimpleGoF1Scores())
					csvCols += "\t" + FRMT.format(bestB2ga.getSimpleGoAnnotationScore());
				if (getSettings().getCalculateAncestryGoF1Scores())
					csvCols += "\t" + FRMT.format(bestB2ga.getAncestryGoAnnotationScore());
				if (getSettings().getCalculateSemSimGoF1Scores())
					csvCols += "\t" + FRMT.format(bestB2ga.getSemSimGoAnnotationScore());
			}
		} else {
			csvCols += "\t\t0\t0.0";
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
				csvCols += "\t";
				if (getSettings().getCalculateSimpleGoF1Scores())
					csvCols += "\t0.0";
				if (getSettings().getCalculateAncestryGoF1Scores())
					csvCols += "\t0.0";
				if (getSettings().getCalculateSemSimGoF1Scores())
					csvCols += "\t0.0";
			}
		}

		return csvCols;
	}

	/**
	 * @param Protein
	 *            prot
	 * @return String - Part of the CSV-Row with the columns Evaluator-Score and
	 *         Reference-Description.
	 */
	public String buildTrainerColumns(Protein prot) {
		// HEADER:
		// \tHRD-Length\tReference-Description\tRef-Lenght\tEvaluation-Score\tDiff-to-bestCompetitor
		String csvCells = "";
		// HRD-Length reference and AHRD's performance:
		if (prot.getEvaluationScoreCalculator().getEvalutionScore() != null) {
			// HRD-Length ref f1score diff-to-best-competitor:
			csvCells += "\t";
			if (prot.getDescriptionScoreCalculator().getHighestScoringBlastResult() != null)
				csvCells += prot.getDescriptionScoreCalculator().getHighestScoringBlastResult().getEvaluationTokens()
						.size();
			else
				csvCells += "0";
			csvCells += "\t" + prot.getEvaluationScoreCalculator().getReferenceDescription().getDescription() + "\t"
					+ prot.getEvaluationScoreCalculator().getReferenceDescription().getTokens().size() + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getEvalutionScore()) + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getEvalScoreMinBestCompScore()) + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getTruePositivesRate()) + "\t"
					+ FRMT.format(prot.getEvaluationScoreCalculator().getFalsePositivesRate());
		} else
			csvCells = "\t\t\t\t\t\t\t";
		return csvCells;
	}

	/**
	 * Append the following columns to the current CSV-Row:
	 * sum_of_all_token_scores, token_high_score, correction_factor, go_score,
	 * lexical_score, rel_bit_score, desc_line_frequency,
	 * max_desc_line_frequency
	 * 
	 * @note: Printing out these values may slow down the running-time, as these
	 *        values are not necessarily stored in memory.
	 * 
	 * @param prot
	 * @return String - Part of the CSV-Row holding the values of the above
	 *         columns.
	 */
	public String buildDescScoreCells(Protein prot) {
		String csvCells = "";
		// Found a high scoring description?
		if (prot.getDescriptionScoreCalculator().getHighestScoringBlastResult() == null) {
			csvCells = "\t\t\t\t\t\t\t\t";
		} else {
			BlastResult hsbr = prot.getDescriptionScoreCalculator().getHighestScoringBlastResult();
			csvCells += "\t" + FRMT.format(prot.getTokenScoreCalculator().sumOfAllTokenScores(hsbr));
			csvCells += "\t" + FRMT.format(prot.getTokenScoreCalculator().getTokenHighScore());
			csvCells += "\t" + FRMT.format(prot.getLexicalScoreCalculator().correctionFactor(hsbr));
			csvCells += "\t" + FRMT.format(prot.getLexicalScoreCalculator().lexicalScore(hsbr));
			csvCells += "\t" + FRMT.format(prot.getDescriptionScoreCalculator().relativeBlastScore(hsbr));
		}
		return csvCells;
	}

	public String buildTokenSetCell(Protein prot) {
		String tokenSetCell = "\t";
		for (String token : prot.getTokenScoreCalculator().getTokenScores().keySet()) {
			tokenSetCell += "[" + token + "->" + FRMT.format(prot.getTokenScoreCalculator().getTokenScores().get(token))
					+ "]";
		}
		return tokenSetCell;
	}

	public String buildBestBlastHitsHeader() {
		String hdr = "";
		for (String blastDb : getSettings().getBlastDatabases()) {
			if (blastDb != null && !blastDb.equals(""))
				hdr += ("\tBest BlastHit against '" + blastDb + "'");
			if (getSettings().isInTrainingMode()) {
				hdr += "\tLength\tEvaluation-Score";
				if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
					hdr += "\tBlast2GO-Annotations '" + blastDb + "'";
					if (getSettings().getCalculateSimpleGoF1Scores())
						hdr += "\tBlast2GO-Annotations-Simple-F-Score '" + blastDb + "'";
					if (getSettings().getCalculateAncestryGoF1Scores())
						hdr += "\tBlast2GO-Annotations-Ancestry-F-Score '" + blastDb + "'";
					if (getSettings().getCalculateSemSimGoF1Scores())
						hdr += "\tBlast2GO-Annotations-SemSim-F-Score '" + blastDb + "'";
				}
			}
		}
		return hdr;
	}

	public String buildBestBlastHitsColumns(Protein prot) {
		String csvRow = "";
		for (String blastDb : getSettings().getBlastDatabases()) {
			if (prot.getEvaluationScoreCalculator().getUnchangedBlastResults().get(blastDb) != null) {
				BlastResult bestBr = prot.getEvaluationScoreCalculator().getUnchangedBlastResults().get(blastDb);
				csvRow += "\t\"" + bestBr.getAccession() + " " + bestBr.getDescription() + "\"";
				if (bestBr.getEvaluationScore() != null) {
					csvRow += "\t" + bestBr.getEvaluationTokens().size() + "\t"
							+ FRMT.format(bestBr.getEvaluationScore());
					if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
						csvRow += "\t" + combineGoTermsToString(bestBr.getGoAnnotations());
						if (getSettings().getCalculateSimpleGoF1Scores())
							csvRow += "\t" + FRMT.format(bestBr.getSimpleGoAnnotationScore());
						if (getSettings().getCalculateAncestryGoF1Scores())
							csvRow += "\t" + FRMT.format(bestBr.getAncestryGoAnnotationScore());
						if (getSettings().getCalculateSemSimGoF1Scores())
							csvRow += "\t" + FRMT.format(bestBr.getSemSimGoAnnotationScore());
					}
				}
			} else {
				csvRow += "\t";
				if (getSettings().isInTrainingMode()) {
					csvRow += "\t0\t0.0";
					if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
						csvRow += "\t";
						if (getSettings().getCalculateSimpleGoF1Scores())
							csvRow += "\t0.0";
						if (getSettings().getCalculateAncestryGoF1Scores())
							csvRow += "\t0.0";
						if (getSettings().getCalculateSemSimGoF1Scores())
							csvRow += "\t0.0";
					}
				}
			}
		}
		return csvRow;
	}

	private String buildReferenceGoAnnotationColumns(Protein prot) {
		String goColumns = "\t"
				+ combineGoTermsToString(prot.getEvaluationScoreCalculator().getReferenceGoAnnoatations());
		if (getSettings().getCalculateSimpleGoF1Scores())
			goColumns += "\t" + FRMT.format(prot.getEvaluationScoreCalculator().getSimpleGoAnnotationScore());
		if (getSettings().getCalculateAncestryGoF1Scores())
			goColumns += "\t" + FRMT.format(prot.getEvaluationScoreCalculator().getAncestryGoAnnotationScore());
		if (getSettings().getCalculateSemSimGoF1Scores())
			goColumns += "\t" + FRMT.format(prot.getEvaluationScoreCalculator().getSemSimGoAnnotationScore());
		return goColumns;
	}

	public String combineGoTermsToString(Set<GOterm> gos) {
		return combineGoTermsToString(gos, ", ");
	}

	public String combineGoTermsToString(Set<GOterm> gos, String seperator) {
		String goLine = "";
		if (gos != null) {
			List<String> sortedGos = new ArrayList<String>();
			for (Iterator<GOterm> goTermIter = gos.iterator(); goTermIter.hasNext();) {
				GOterm term = goTermIter.next();
				sortedGos.add(term.getAccession());
			}
			Collections.sort(sortedGos);
			for (Iterator<String> iter = sortedGos.iterator(); iter.hasNext();) {
				String term = iter.next();
				goLine += term;
				if (iter.hasNext())
					goLine += seperator;
			}

		}
		return goLine;
	}

}
