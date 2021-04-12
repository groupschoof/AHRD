package ahrd.view;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import ahrd.model.BlastResult;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;

public abstract class OutputWriter {

	/**
	 * Format decimal numbers to three digits after decimal-point and leading
	 * zero, if number is smaller than zero.
	 */
	public static final DecimalFormat FRMT = new DecimalFormat("#,###0.###");
	static {
        final DecimalFormatSymbols dcs = FRMT.getDecimalFormatSymbols();
        dcs.setDecimalSeparator('.');
        dcs.setGroupingSeparator(',');
        dcs.setNaN("NaN");
        FRMT.setDecimalFormatSymbols(dcs);
	}
	
	private Collection<Protein> proteins;

	public OutputWriter(Collection<Protein> proteins) {
		setProteins(proteins);
	}

	public abstract void writeOutput() throws IOException;

	public String buildDescriptionLine(Protein protein, String separator) {
		String descLine = protein.getAccession() + separator;
		// Blast-Results
		if (protein.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult() != null) {
			BlastResult br = protein.getDescriptionScoreCalculator()
					.getHighestScoringBlastResult();
			descLine += br.getAccession() + separator + qualityCode(protein)
					+ separator + br.getDescription() + separator;
		} else {
			// Maintain Table's Column-Structure, if writing tab delimited
			// values:
			if (separator.equals("\t"))
				descLine += "\t\tUnknown protein\t";
			else
				descLine += "Unknown protein";
		}
		descLine += separator;
		// Gene-Ontology-Results:
		descLine += combineGoTermStrings(protein.getGoResults());
		
		return descLine;
	}
	
	public String combineGoTermStrings(Set<String> gos) {
		return combineGoTermStrings(gos, ", ");
	}
	
	public String combineGoTermStrings(Set<String> gos, String separator) {
		String goLine = "";
		List<String> sortedGOs = new ArrayList<String>(gos);
		Collections.sort(sortedGOs);
		for (Iterator<String> i = sortedGOs.iterator(); i.hasNext();) {
			String gor = i.next();
			goLine += gor;
			if (i.hasNext())
				goLine += separator;
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
