package ahrd.view;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import ahrd.model.BlastResult;
import ahrd.model.InterproResult;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;

public abstract class AbstractOutputWriter implements OutputWriter {

	/**
	 * Format decimal numbers to three digits after decimal-point and leading
	 * zero, if number is smaller than zero.
	 */
	public static final DecimalFormat FRMT = new DecimalFormat("#,###0.###");

	private Collection<Protein> proteins;

	public AbstractOutputWriter(Collection<Protein> proteins) {
		FRMT.getDecimalFormatSymbols().setDecimalSeparator('.');
		FRMT.getDecimalFormatSymbols().setGroupingSeparator(',');
		setProteins(proteins);
	}

	public abstract void writeOutput() throws IOException;

	public String buildDescriptionLine(Protein protein, String seperator) {
		String descLine = protein.getAccession() + seperator;
		// Blast-Results
		if (protein.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult() != null) {
			BlastResult br = protein.getDescriptionScoreCalculator()
					.getHighestScoringBlastResult();
			descLine += br.getAccession() + seperator + qualityCode(protein)
					+ seperator + br.getDescription() + seperator;
		} else {
			// Maintain Table's Column-Structure, if writing tab delimited
			// values:
			if (seperator.equals("\t"))
				descLine += "\t\tUnknown protein\t";
			else
				descLine += "Unknown protein";
		}
		// Interpro
		List<InterproResult> sortedIprs = new ArrayList<InterproResult>(
				protein.getInterproResults());
		Collections.sort(sortedIprs);
		for (Iterator<InterproResult> i = sortedIprs.iterator(); i.hasNext();) {
			InterproResult ipr = i.next();
			descLine += ipr.getId() + " (" + ipr.getName() + ")";
			if (i.hasNext())
				descLine += ", ";
		}
		descLine += seperator;
		// Gene-Ontology-Results:
		descLine += combineGoTermStrings(protein.getGoResults());
		
		return descLine;
	}
	
	public String combineGoTermStrings(Set<String> gos) {
		return combineGoTermStrings(gos, ", ");
	}
	
	public String combineGoTermStrings(Set<String> gos, String seperator) {
		String goLine = "";
		List<String> sortedGOs = new ArrayList<String>(gos);
		Collections.sort(sortedGOs);
		for (Iterator<String> i = sortedGOs.iterator(); i.hasNext();) {
			String gor = i.next();
			goLine += gor;
			if (i.hasNext())
				goLine += seperator;
		}
		return goLine;
	}
	
	/**
	 * Four Positions. Each gets a '*', if...
	 * 
	 * Position One: BitScore > 50 and EValue < 0.1
	 * 
	 * Position Two: Overlap > 60%
	 * 
	 * Position Three: DescriptionScore >= 0.5
	 * 
	 * Position Four: DescriptionLine shares tokens with predicted
	 * Gene-Ontology-Terms
	 * 
	 * @return the quality code
	 */
	public String qualityCode(Protein p) {
		BlastResult hsbr = p.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult();
		String qc = "";
		// Position 1
		qc += (hsbr.getBitScore() > 50.0 && hsbr.getEValue() < 0.1) ? "*" : "-";
		// Position 2
		qc += (TokenScoreCalculator.overlapScore(hsbr.getQueryStart(),
				hsbr.getQueryEnd(), p.getSequenceLength(),
				hsbr.getSubjectStart(), hsbr.getSubjectEnd(),
				hsbr.getSubjectLength()) > 0.6) ? "*" : "-";
		// Position 3
		qc += (p.getDescriptionScoreCalculator().getDescriptionHighScore() >= 0.5) ? "*"
				: "-";
		// Internal DescriptionScore:
		// qc += "[" + FRMT.format(hsbr.getDescriptionScore()) + "]";

		return qc;
	}

	public Collection<Protein> getProteins() {
		return proteins;
	}

	public void setProteins(Collection<Protein> proteins) {
		this.proteins = proteins;
	}

}
