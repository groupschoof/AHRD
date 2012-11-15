package ahrd.view;

import static ahrd.controller.Settings.getSettings;
import static ahrd.model.TokenScoreCalculator.overlapScore;
import static ahrd.view.AbstractOutputWriter.FRMT;
import static ahrd.view.AbstractOutputWriter.geneOntologyAnnotations;
import static ahrd.view.AbstractOutputWriter.interproResult;
import static ahrd.view.AbstractOutputWriter.qualityCode;

import java.util.HashMap;
import java.util.Map;

import ahrd.model.BlastResult;
import ahrd.model.Protein;

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
	public String HRDLength;
	public String ReferenceDescription;
	public String RefLength;
	public String EvaluationScore;
	public String DifftobestCompetitor;
	public String TPR;
	public String FPR;
	public Map<String, BlastResultOutput> Best_BlastHits = new HashMap<String, BlastResultOutput>();
	public String Sum_TokenScores;
	public String TokenHighScore;
	public String CorrectionFactor;
	public String GOScore;
	public String OverlapScore;
	public String LexicalScore;
	public String RelativeBitScore;
	public String DescriptionLineFrequency;
	public String Max_DescLineFreq;
	public String PatternFactor;
	public String vectorSpaceModel;
	public String ProteinDomainWeightVector;
	public String HRDDomainArchitectureDomainWeightVector;
	public String HRDDomainArchitectureSimilarityScore;

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
		BlastResult br = prot.getDescriptionScoreCalculator()
				.getHighestScoringBlastResult();
		this.ProteinAccession = prot.getAccession();
		if (br != null) {
			this.BlastHitAccession = br.getAccession();
			this.AHRDQualityCode = qualityCode(prot);
			this.HumanReadableDescription = br.getDescription();
		}
		if (prot.getInterproResults() != null
				&& !prot.getInterproResults().isEmpty()) {
			this.InterproID = interproResult(prot);
		} else
			this.InterproID = NA;
		if (prot.getGoResults() != null && !prot.getGoResults().isEmpty()) {
			this.GeneOntologyID = geneOntologyAnnotations(prot);
		}
		if (getSettings().isInTrainingMode()) {
			// HRD-Length\tReference-Description\tRef-Length\tEvaluation-Score\tDiff-to-bestCompetitor\tTPR\tFPR"
			if (br != null && br.getTokens() != null)
				this.HRDLength = format(br.getTokens().size());
			this.ReferenceDescription = prot.getEvaluationScoreCalculator()
					.getReferenceDescription().getDescription();
			this.RefLength = format(prot.getEvaluationScoreCalculator()
					.getReferenceDescription().getTokens().size());
			this.EvaluationScore = format(prot.getEvaluationScoreCalculator()
					.getEvalutionScore());
			this.DifftobestCompetitor = format(prot
					.getEvaluationScoreCalculator()
					.getEvalScoreMinBestCompScore());
			this.TPR = format(prot.getEvaluationScoreCalculator()
					.getTruePositivesRate());
			this.FPR = format(prot.getEvaluationScoreCalculator()
					.getFalsePositivesRate());
		}
		if (getSettings().isWriteDomainArchitectureSimilarityScoresToOutput()) {
			this.ProteinDomainWeightVector = saveToString(prot
					.getDomainWeights());
			if (br != null)
				this.HRDDomainArchitectureSimilarityScore = format(br
						.getDomainSimilarityScore());
			else
				this.HRDDomainArchitectureSimilarityScore = NA;
			this.vectorSpaceModel = saveToString(prot
					.getDomainScoreCalculator().getVectorSpaceModel());
			if (br != null)
				this.HRDDomainArchitectureDomainWeightVector = saveToString(br
						.getDomainWeights());
			else
				this.HRDDomainArchitectureDomainWeightVector = NA;
		}
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			for (String blastDb : getSettings().getBlastDatabases()) {
				BlastResult iterBr = prot.getEvaluationScoreCalculator()
						.getUnchangedBlastResults().get(blastDb);
				this.Best_BlastHits.put(blastDb, new BlastResultOutput(iterBr));
			}
		}
		if (getSettings().doWriteHRDScoresToOutput()) {
			this.Sum_TokenScores = br != null ? format(prot
					.getTokenScoreCalculator().sumOfAllTokenScores(br)) : NA;
			this.TokenHighScore = br != null ? format(prot
					.getTokenScoreCalculator().getTokenHighScore()) : NA;
			this.CorrectionFactor = br != null ? format(prot
					.getLexicalScoreCalculator().correctionFactor(br)) : NA;
			this.GOScore = br != null ? format(prot.getLexicalScoreCalculator()
					.geneOntologyScore(br)) : NA;
			this.OverlapScore = br != null ? format(overlapScore(
					br.getQueryStart(), br.getQueryEnd(),
					prot.getSequenceLength(), br.getSubjectStart(),
					br.getSubjectEnd(), br.getSubjectLength())) : NA;
			this.LexicalScore = br != null ? format(prot
					.getLexicalScoreCalculator().lexicalScore(br)) : NA;
			this.RelativeBitScore = br != null ? format(prot
					.getDescriptionScoreCalculator().relativeBlastScore(br))
					: NA;
			this.DescriptionLineFrequency = br != null ? format(prot
					.getDescriptionScoreCalculator()
					.getDescLinePatternFrequencies().get(br.patternize())) : NA;
			this.Max_DescLineFreq = br != null ? format(prot
					.getDescriptionScoreCalculator()
					.getMaxDescriptionLineFrequency()) : NA;
			this.PatternFactor = br != null ? format(prot
					.getDescriptionScoreCalculator().patternFactor(br)) : NA;
		}
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
