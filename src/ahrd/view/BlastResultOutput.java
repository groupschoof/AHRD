package ahrd.view;

import static ahrd.controller.Settings.getSettings;
import static ahrd.view.ProteinOutput.saveToString;
import static ahrd.view.ProteinOutput.format;
import java.util.Set;

import ahrd.model.BlastResult;

/**
 * This class is meant as a mediator between output writers and the model
 * classes. It extracts and presents those values of a query protein that have
 * to go into the output. It enables usage of libraries that convert a given
 * Java Object to XML or JSON.
 * 
 * @Note: This class does not use getters and setters, but declares public
 *        fields to ease access from other libraries.
 * 
 * @author Asis Hallab, Kathrin Klee, Mythri Bangalore
 */
public class BlastResultOutput {

	public String accession = ProteinOutput.NA;
	public String description = ProteinOutput.NA;
	public String length;
	public String domainWeightsVector;
	public String domainSimilarityScore;

	public BlastResultOutput(BlastResult blastResult) {
		if (blastResult != null) {
			this.accession = blastResult.getAccession();
			this.description = blastResult.getDescription();
		}
		if (getSettings().isInTrainingMode()) {
			if (blastResult != null)
				setLength(blastResult.getTokens());
			else
				this.length = ProteinOutput.NA;
		}
		if (getSettings().isWriteDomainArchitectureSimilarityScoresToOutput()) {
			if (blastResult != null) {
				this.domainWeightsVector = saveToString(blastResult
						.getDomainWeights());
				this.domainSimilarityScore = format(blastResult
						.getDomainSimilarityScore());
			} else {
				this.domainWeightsVector = ProteinOutput.NA;
				this.domainSimilarityScore = ProteinOutput.NA;
			}
		}
	}

	public void setLength(Set<String> tokens) {
		if (tokens != null && tokens.size() > 0)
			this.length = new Integer(tokens.size()).toString();
		else
			this.length = "0";
	}
}
