package ahrd.view;

import static ahrd.controller.Settings.getSettings;
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

	public BlastResultOutput(BlastResult blastResult) {
		if (blastResult != null) {
			this.accession = blastResult.getAccession();
			this.description = blastResult.getDescription();
		}
		if (getSettings().isInTrainingMode())
			setLength(blastResult);
	}

	public void setLength(BlastResult blastResult) {
		if (blastResult != null && blastResult.getTokens() != null) {
			if (blastResult.getTokens() != null
					&& blastResult.getTokens().size() > 0)
				this.length = new Integer(blastResult.getTokens().size())
						.toString();
		} else
			this.length = "0";
	}
}
