package ahrd.view;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import ahrd.model.BlastResult;
import ahrd.model.GeneOntologyResult;
import ahrd.model.InterproResult;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;

public abstract class AbstractOutputWriter implements IOutputWriter {

	/**
	 * Format decimal numbers to three digits after decimal-point and leading
	 * zero, if number is smaller than zero.
	 */
	public static final DecimalFormat FRMT = new DecimalFormat("#,###0.###");

	private Collection<Protein> proteins;

	public AbstractOutputWriter(Collection<Protein> proteins) {
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
		descLine += interproResult(protein);
		descLine += seperator;
		// Gene-Ontology-Results:
		descLine += geneOntologyAnnotations(protein);
		return descLine;
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
	public static String qualityCode(Protein p) {
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
		// Position 4
		qc += (p.getLexicalScoreCalculator().geneOntologyScore(hsbr) > 0.0) ? "*"
				: "-";
		// Internal DescriptionScore:
		qc += "[" + FRMT.format(hsbr.getDescriptionScore()) + "]";
		
		return qc;
	}

	public Collection<Protein> getProteins() {
		return proteins;
	}

	public void setProteins(Collection<Protein> proteins) {
		this.proteins = proteins;
	}

	public static String interproResult(Protein protein) {
		String descLine = "";
		List<InterproResult> sortedIprs = new ArrayList<InterproResult>(
				protein.getInterproResults());
		Collections.sort(sortedIprs);
		for (Iterator<InterproResult> i = sortedIprs.iterator(); i.hasNext();) {
			InterproResult ipr = i.next();
			descLine += ipr.getId() + " (" + ipr.getName() + ")";
			if (i.hasNext())
				descLine += ", ";
		}
		return descLine;
	}

	public static String geneOntologyAnnotations(Protein protein) {
		String descLine = "";
		List<GeneOntologyResult> sortedGOs = new ArrayList<GeneOntologyResult>(
				protein.getGoResults());
		Collections.sort(sortedGOs);
		for (Iterator<GeneOntologyResult> i = sortedGOs.iterator(); i.hasNext();) {
			GeneOntologyResult gor = i.next();
			descLine += gor.getAcc() + " (" + gor.getName() + ")";
			if (i.hasNext())
				descLine += ", ";
		}
		return descLine;
	}

}
