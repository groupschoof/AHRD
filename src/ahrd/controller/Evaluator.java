package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.IOException;

import ahrd.exception.MissingAccessionException;
import ahrd.model.Blast2GoAnnot;
import ahrd.model.Protein;
import ahrd.model.ReferenceDescription;
import ahrd.view.IOutputWriter;

public class Evaluator extends AHRD {

	public Evaluator(String pathToInputYml) throws IOException {
		super(pathToInputYml);
	}

	public void setupReferences() throws IOException, MissingAccessionException {
		String[] fastaEntries = getSettings().getReferencesFasta().split(">");
		for (String fastaEntry : fastaEntries) {
			if (fastaEntry != null && !fastaEntry.trim().equals("")) {
				ReferenceDescription rd = ReferenceDescription
						.constructFromFastaEntry(fastaEntry.trim());
				Protein p = getProteins().get(rd.getAccession());
				if (p == null)
					throw new MissingAccessionException(
							"Could not find Protein for Accession '"
									+ rd.getAccession() + "'");
				p.getEvaluationScoreCalculator().setReferenceDescription(rd);
			}
		}
	}

	public void setupBlast2GoAnnots() throws IOException,
			MissingAccessionException {
		if (getSettings().getPathToBlast2GoAnnotations() != null
				&& !getSettings().getPathToBlast2GoAnnotations().equals("")) {
			for (String blast2GoResultEntry : getSettings()
					.getBlast2GoAnnotations()) {
				Blast2GoAnnot b2ga = Blast2GoAnnot
						.fromBlast2GoEntry(blast2GoResultEntry);
				Protein p = getProteins().get(b2ga.getAccession());
				if (p == null)
					throw new MissingAccessionException(
							"Could not find Protein for Accession '"
									+ b2ga.getAccession() + "'");
				p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2ga);
			}
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out
				.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.Evaluator input.yml\n");

		try {
			Evaluator evaluator = new Evaluator(args[0]);
			evaluator.setup(false); // false -> Don't log memory and time-usages
			evaluator.setupReferences();
			// Blast2GO is another competitor in the field of annotation of
			// predicted Proteins. AHRD might be compared with B2Gs performance:
			evaluator.setupBlast2GoAnnots();
			// Iterate over all Proteins and assign the best scoring Human
			// Readable Description
			evaluator.assignHumanReadableDescriptions();
			// Evaluate AHRD's performance for each Protein:
			evaluator.calculateEvaluationScores();
			// If requested, calculate the highest possibly achievable
			// evaluation score:
			if (getSettings().doFindHighestPossibleEvaluationScore())
				evaluator.findHighestPossibleEvaluationScores();
			// Generate Output:
			IOutputWriter ow = initializeOutputWriter(evaluator.getProteins()
					.values());
			ow.writeOutput();
			System.out.println("Written output into:\n"
					+ getSettings().getPathToOutput());
		} catch (Exception e) {
			System.err.println("We are sorry, an unexpected ERROR occurred:");
			e.printStackTrace(System.err);
		}

	}

	/**
	 * Calculates the evaluation-score for the description assigned by AHRD and
	 * all found descriptions coming from competitive methods. As for now those
	 * methods are the best Blast-Hits of each searched Blast-Database. @NOTE:
	 * Best Blast-Hits have the highest Bit-Scores. <br />
	 * After having evaluated the objective function for each Protein,
	 * calculates the average-evaluation-score and sets it in the current
	 * Settings.
	 */
	public void calculateEvaluationScores() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator().assignEvlScrsToCompetitors();
		}
	}

	/**
	 * Calculates the evaluation-score as defined in function
	 * calculateEvaluationScores() for each Protein's Blast-Results'
	 * Description-Lines and stores the maximum possible score each Protein's
	 * EvaluationScoreCalculator. This is used to get more accurate information
	 * of how well AHRD performs.
	 */
	public void findHighestPossibleEvaluationScores() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator()
					.findHighestPossibleEvaluationScore();
		}
	}

}
