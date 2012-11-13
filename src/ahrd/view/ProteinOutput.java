package ahrd.view;

import java.util.HashMap;
import java.util.Map;

import ahrd.model.BlastResult;
import ahrd.model.Protein;
import static ahrd.view.AbstractOutputWriter.qualityCode;
import static ahrd.view.AbstractOutputWriter.interproResult;
import static ahrd.view.AbstractOutputWriter.geneOntologyAnnotations;
import static ahrd.controller.Settings.getSettings;
import static ahrd.view.AbstractOutputWriter.FRMT;

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
public class ProteinOutput {

	/**
	 * The following fields are named as closely as possible as their
	 * corresponding names in the output:
	 */
	public String ProteinAccession;
	public String BlastHitAccession = NA;
	public String AHRDQualityCode = NA;
	public String HumanReadableDescription = "unknown protein";
	public String InterproID = NA;
	public String GeneOntologyID = NA;
	public String HRDLength = NA;
	public String ReferenceDescription;
	public String RefLength;
	public String EvaluationScore;
	public String DifftobestCompetitor;
	public String TPR;
	public String FPR;
	public String ProteinDomainWeightVector;
	public String HRDDomainArchitectureSimilarityScore;
	public Map<String, BlastResultOutput> Best_BlastHits = new HashMap<String, BlastResultOutput>();
	public String Length;
	public String BHDomainWeightVector;
	public String BHDomainArchitectureSimilarityScore;
	public String Sum_TokenScores;
	public String TokenHighScore;
	public String CorrectionFactor;
	public String GOScore;
	public String LexicalScore;
	public String RelativeBitScore;
	public String DescriptionLineFrequency;
	public String Max_DescLineFreq;
	public String PatternFactor;
	public String vectorSpaceModel;
	public String HRDDomainArchitectureDomainWeightVector;

	/**
	 * NULL will be translated into 'NA'
	 */
	public static final String NA = "NA";

	/**
	 * Empty constructor.
	 */
	public ProteinOutput() {
		super();
	}

	/**
	 * @param prot
	 */
	public ProteinOutput(Protein prot) {
		super();
		System.out.println("1");
		BlastResult br = prot.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult();
		System.out.println("2");
		this.ProteinAccession = prot.getAccession();
		System.out.println("3");
		if (br != null) {
			this.BlastHitAccession = br.getAccession();
			System.out.println("4");
			this.AHRDQualityCode = qualityCode(prot);
			System.out.println("5");
			this.HumanReadableDescription = br.getDescription();
			System.out.println("6");
		}
		if (prot.getInterproResults() != null
				&& !prot.getInterproResults().isEmpty()) {
			this.InterproID = interproResult(prot);
			System.out.println("7");
		} else
			this.InterproID = NA;
		if (prot.getGoResults() != null && !prot.getGoResults().isEmpty()) {
			this.GeneOntologyID = geneOntologyAnnotations(prot);
			System.out.println("8");
		}
		if (getSettings().isInTrainingMode()) {
			// HRD-Length\tReference-Description\tRef-Length\tEvaluation-Score\tDiff-to-bestCompetitor\tTPR\tFPR"
			if (br != null && br.getTokens() != null)
				this.HRDLength = format(br.getTokens().size());
			System.out.println("9");
			this.ReferenceDescription = prot.getEvaluationScoreCalculator()
					.getReferenceDescription().getDescription();
			System.out.println("10");
			this.RefLength = format(prot.getEvaluationScoreCalculator()
					.getReferenceDescription().getTokens().size());
			System.out.println("11");
			this.EvaluationScore = format(prot.getEvaluationScoreCalculator()
					.getEvalutionScore());
			System.out.println("12");
			this.DifftobestCompetitor = format(prot
					.getEvaluationScoreCalculator()
					.getEvalScoreMinBestCompScore());
			System.out.println("13");
			this.TPR = format(prot.getEvaluationScoreCalculator()
					.getTruePositivesRate());
			System.out.println("14");
			this.FPR = format(prot.getEvaluationScoreCalculator()
					.getFalsePositivesRate());
			System.out.println("15");
		}
		System.out.println("16_A");
		if (getSettings().isWriteDomainArchitectureSimilarityScoresToOutput()) {
			this.ProteinDomainWeightVector = saveToString(prot
					.getDomainWeights());
			System.out.println("17");
			if (br != null)
				this.HRDDomainArchitectureSimilarityScore = format(br
						.getDomainSimilarityScore());
			else
				this.HRDDomainArchitectureSimilarityScore = NA;
			System.out.println("18");
			this.vectorSpaceModel = saveToString(prot
					.getDomainScoreCalculator().getVectorSpaceModel());
			System.out.println("19");
			if (br != null)
				this.HRDDomainArchitectureDomainWeightVector = saveToString(br
						.getDomainWeights());
			else
				this.HRDDomainArchitectureDomainWeightVector = NA;
			System.out.println("20");
		}
		System.out.println("20_A");
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			for (String blastDb : getSettings().getBlastDatabases()) {
				BlastResult iterBr = prot.getEvaluationScoreCalculator()
						.getUnchangedBlastResults().get(blastDb);
				this.Best_BlastHits.put(blastDb, new BlastResultOutput(iterBr));
				System.out.println("21");
			}
		}
		System.out.println("22");
		if (getSettings().doWriteHRDScoresToOutput() && br != null) {
			this.Sum_TokenScores = format(prot.getTokenScoreCalculator()
					.sumOfAllTokenScores(br));
			System.out.println("23");
			this.TokenHighScore = format(prot.getTokenScoreCalculator()
					.getTokenHighScore());
			System.out.println("24");
			this.CorrectionFactor = format(prot.getLexicalScoreCalculator()
					.correctionFactor(br));
			System.out.println("25");
			this.GOScore = format(prot.getLexicalScoreCalculator()
					.geneOntologyScore(br));
			System.out.println("26");
			this.LexicalScore = format(prot.getLexicalScoreCalculator()
					.lexicalScore(br));
			System.out.println("27");
			this.RelativeBitScore = format(prot.getDescriptionScoreCalculator()
					.relativeBlastScore(br));
			System.out.println("28");
			this.DescriptionLineFrequency = format(prot
					.getDescriptionScoreCalculator()
					.getDescLinePatternFrequencies().get(br.patternize()));
			System.out.println("29");
			this.Max_DescLineFreq = format(prot.getDescriptionScoreCalculator()
					.getMaxDescriptionLineFrequency());
			System.out.println("30");
			this.PatternFactor = format(prot.getDescriptionScoreCalculator()
					.patternFactor(br));
			System.out.println("31");
		}
		System.out.println("32");
	}

	public static String saveToString(Object o) {
		return o != null ? o.toString() : NA;
	}

	public static String saveToString(int o) {
		return new Integer(o).toString();
	}

	public static String format(Number n) {
		return n != null ? FRMT.format(n) : NA;
	}
}
