package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Vector;

import ahrd.exception.MissingInterproResultException;

/**
 * We estimate protein similarity based on the formulae given in the article
 * "Protein comparison at the domain architecture level" by Lee and Lee
 * (http://www.biomedcentral.com/1471-2105/10/S15/S5). The more similar a
 * described protein is to the query the more appropriate its description. Hence
 * we introduce a new Domain-Score on the level of BlastResults to reflect this
 * concept.
 * 
 * @author hallab, bangalore, klee
 */
public class DomainScoreCalculator {

	/**
	 * Result from parsing SIMAP's concatenated and preprocessed feature-files.
	 * Preprocessing involves substitution of SIMAP-Hashes with the original
	 * Protein-Accessions.
	 * 
	 * @note: See awk-scripts in directory helper_scripts.
	 */
	private static Map<String, Set<String>> blastResultAccessionsToInterproIds = new HashMap<String, Set<String>>();

	/**
	 * To enable calculation of domain-architecture scores, we need to know the
	 * concrete architecture of proteins of significant similarity (BLAST
	 * results). In this memory database we store the Pfam domains.
	 */
	private static Map<String, Set<String>> blastResultAccessionsToPfamIds = new HashMap<String, Set<String>>();

	private Protein protein;
	private SortedSet<String> vectorSpaceModel;
	private Map<String, Double> cumulativeTokenDomainSimilarityScores = new HashMap<String, Double>();
	private Double totalTokenDomainSimilarityScore = 0.0;

	/**
	 * 
	 * 1.) Construct the vector space model for the Protein and its BlastResults
	 * 2.) Construct the domain-weights vector for the Protein itself 3.) ...and
	 * all of its BlastResults
	 * 
	 * @param prot
	 * @throws MissingInterproResultException
	 */
	public static void constructDomainWeightVectors(Protein prot)
			throws MissingInterproResultException {

		// Vector Space Model of all distinct annotated Interpro-Entities:
		SortedSet<String> vsm = constructVectorSpaceModel(prot);

		// Domain-Weight Vector for the Protein itself:
		List<Double> prVec = new Vector<Double>();
		for (Iterator<String> it = vsm.iterator(); it.hasNext();) {
			String domainAccession = it.next();
			prVec.add(getDomainWeight(prot, domainAccession));
		}
		// Set the results:
		prot.setDomainWeights(prVec);

		// Domain-Weight Vector for all Protein's BlastResults:
		for (String blastDb : prot.getBlastResults().keySet()) {
			for (BlastResult br : prot.getBlastResults().get(blastDb)) {
				List<Double> brVec = new Vector<Double>();
				for (String domainAccession : vsm) {
					brVec.add(getDomainWeight(br, domainAccession));
				}
				// Set the results:
				br.setDomainWeights(brVec);
			}
		}
	}

	/**
	 * Calculates the cosine of angle between the two argument vectors as a
	 * measure of their similarity: sim(x,y) = dot-product(x,y) / (||x||*||y||).
	 * 
	 * For any or both vectors equaling the origin, this function returns not
	 * NaN, but zero.
	 * 
	 * @param x
	 * @param y
	 * @return sim(x,y)
	 */
	public static Double domainWeightSimilarity(List<Double> prVec,
			List<Double> brVec) {
		Double dotProduct = 0.0;
		for (int i = 0; i < prVec.size(); i++) {
			dotProduct += (prVec.get(i) * brVec.get(i));
		}

		Double magPr = 0.0;
		Double magBr = 0.0;
		for (int i = 0; i < prVec.size(); i++) {
			magPr += Math.pow(prVec.get(i), 2);
			magBr += Math.pow(brVec.get(i), 2);
		}
		Double magnitude = Math.sqrt(magPr) * Math.sqrt(magBr);
		Double dws = dotProduct / magnitude;

		if (dws.equals(Double.NaN))
			dws = 0.0;

		return dws;
	}

	/**
	 * Gathers all distinct InterproIDs assigned to the Protein and its
	 * BlastResults, than constructs a sorted Set of them to be used as a
	 * definition for the vector space model.
	 * 
	 * @Note: Depending on AHRD's input, it either bases the vector space model
	 *        on InterPro or Pfam annotations.
	 * 
	 * @param Protein
	 *            prot
	 * @return SortedSet<String> of the respective InterPro or Pfam identifiers
	 *         in their natural order
	 */
	public static SortedSet<String> constructVectorSpaceModel(Protein prot) {
		SortedSet<String> vectorSpaceModel = new TreeSet<String>();
		// Add all Protein's domain annotation:
		vectorSpaceModel.addAll(getDomainAnnotation(prot));
		// Add all domain annotation of all Protein's BlastResults:
		for (String blastDb : prot.getBlastResults().keySet()) {
			for (BlastResult br : prot.getBlastResults().get(blastDb)) {
				vectorSpaceModel.addAll(getDomainAnnotation(br));
			}
		}
		return vectorSpaceModel;
	}

	/**
	 * Extracts either the Pfam or the InterPro annotations from the argument,
	 * depending which is set to be used in AHRD's input.
	 * 
	 * @param Protein
	 *            prot
	 * @return Set<String> domain annotation
	 */
	public static Set<String> getDomainAnnotation(Protein prot) {
		Set<String> domainAnnotation = new HashSet<String>();
		if (getSettings()
				.isDomainArchitectureSimilarityBasedOnPfamAnnotations())
			domainAnnotation = prot.getPfamResults();
		else
			for (InterproResult ipr : prot.getInterproResults()) {
				domainAnnotation.add(ipr.getId());
			}
		return domainAnnotation;
	}

	/**
	 * Extracts either the Pfam or the InterPro annotations from the argument,
	 * depending which is set to be used in AHRD's input.
	 * 
	 * @param Protein
	 *            prot
	 * @return Set<String> domain annotation, <em>empty</em> if none found.
	 */
	public static Set<String> getDomainAnnotation(BlastResult blastResult) {
		Set<String> domainAnnotation = new HashSet<String>();
		Set<String> blastResultAnnotation = null;
		if (getSettings()
				.isDomainArchitectureSimilarityBasedOnPfamAnnotations())
			blastResultAnnotation = getBlastResultAccessionsToPfamIds().get(
					blastResult.getAccession());
		else
			blastResultAnnotation = getBlastResultAccessionsToInterproIds()
					.get(blastResult.getAccession());
		// Validate:
		if (blastResultAnnotation != null && !blastResultAnnotation.isEmpty())
			domainAnnotation = blastResultAnnotation;

		return domainAnnotation;
	}

	/**
	 * Reads out the domain weight appropriate for the argument. This is either
	 * the domain weight recorded for the conserved protein domain or zero if
	 * this domain has not been annotated for the argument protein.
	 * 
	 * @TODO: Refactor to Interface with method getDomainWeight implemented by
	 *        both classes BlastResult and Protein.
	 * 
	 * @param br
	 * @param domainAccession
	 * @return Double - 0.0 is returned to avoid NullPointerExceptions, in case
	 *         the InterPro Entry could be found, but it had a NULL domain
	 *         weight.
	 * @throws MissingInterproResultException
	 */
	public static Double getDomainWeight(Protein prot, String domainAccession)
			throws MissingInterproResultException {
		Double dw = 0.0;
		if (getSettings()
				.isDomainArchitectureSimilarityBasedOnPfamAnnotations()) {
			if (prot.getPfamResults().contains(domainAccession)) {
				dw = InterproResult.getPfamDomainWeights().get(domainAccession);
				if (dw == null)
					throw new MissingInterproResultException(
							"Could not find domain weight for Pfam Entry '"
									+ domainAccession + "'in memory database.");
			}
		} else {
			InterproResult ipr = InterproResult.getInterproDb().get(
					domainAccession);
			if (ipr == null)
				throw new MissingInterproResultException(
						"Could not find Interpro-Entry '" + domainAccession
								+ "' in memory database.");
			if (prot.getInterproResults().contains(ipr))
				dw = ipr.getDomainWeight();
		}
		return dw;
	}

	/**
	 * Reads out the domain weight appropriate for the argument. This is either
	 * the domain weight recorded for the conserved protein domain or zero if
	 * this domain has not been annotated for the argument protein.
	 * 
	 * @TODO: Refactor to Interface with method getDomainWeight implemented by
	 *        both classes BlastResult and Protein.
	 * 
	 * @param br
	 * @param domainAccession
	 * @return Double - 0.0 is returned to avoid NullPointerExceptions, in case
	 *         the InterPro Entry could be found, but it had a NULL domain
	 *         weight.
	 * @throws MissingInterproResultException
	 */
	public static Double getDomainWeight(BlastResult br, String domainAccession)
			throws MissingInterproResultException {
		Double dw = 0.0;
		if (getSettings()
				.isDomainArchitectureSimilarityBasedOnPfamAnnotations()) {
			if (getBlastResultAccessionsToPfamIds().containsKey(
					br.getAccession())
					&& getBlastResultAccessionsToPfamIds().get(
							br.getAccession()).contains(domainAccession)) {
				dw = InterproResult.getPfamDomainWeights().get(domainAccession);
				System.out.println("BlastResult: " + br.getAccession()
						+ ", DomainWeights-Db: "
						+ InterproResult.getPfamDomainWeights());
				if (dw == null)
					throw new MissingInterproResultException(
							"Could not find domain weight for Pfam Entry '"
									+ domainAccession + "'in memory database.");
			}
		} else {
			InterproResult ipr = InterproResult.getInterproDb().get(
					domainAccession);
			if (ipr == null)
				throw new MissingInterproResultException(
						"Could not find Interpro-Entry '" + domainAccession
								+ "' in memory database.");
			if (getBlastResultAccessionsToInterproIds().containsKey(
					br.getAccession())
					&& getBlastResultAccessionsToInterproIds().get(
							br.getAccession()).contains(domainAccession))
				dw = ipr.getDomainWeight();
		}
		return dw;
	}

	public DomainScoreCalculator(Protein protein) {
		super();
		setProtein(protein);
	}

	/**
	 * After initialization of BlastResults Interpro-Annotations and having
	 * selected all BlastResults to be description candidates, now the
	 * DomainSimilarityScores can be computed.
	 * <ul>
	 * <li>Generate the vector space model (VSM)</li>
	 * <li>For the query protein itself and each BlastResult (Description
	 * Candidate) compute its Domain Weights Vector</li>
	 * <li>Finally for each BlastResult compute its DomainSimilarityScore as
	 * sim(query-protein, blast-result)</li>
	 * <li>Trigger computation of cumulative token domain similarity scores and
	 * total token domain similarity score in the query protein's
	 * TokenScoreCalculator</li>
	 * <li>Trigger memorization off the maximum found domain similarity score</li>
	 * </ul>
	 * 
	 * @param Protein
	 *            prot
	 * @throws MissingInterproResultException
	 */
	public void computeDomainSimilarityScores()
			throws MissingInterproResultException {
		setVectorSpaceModel(constructVectorSpaceModel(getProtein()));
		constructDomainWeightVectors(getProtein());
		for (String blastDb : getProtein().getBlastResults().keySet()) {
			for (BlastResult br : getProtein().getBlastResults().get(blastDb)) {
				br.setDomainSimilarityScore(domainWeightSimilarity(getProtein()
						.getDomainWeights(), br.getDomainWeights()));
				getProtein().getTokenScoreCalculator()
						.measureCumulativeDomainSimilarityScores(br);
				getProtein().getTokenScoreCalculator()
						.measureTotalDomainSimilarityScore(br);
				getProtein().getDescriptionScoreCalculator()
						.measureMaxDomainSimilarityScore(br);
			}
		}
	}

	public Protein getProtein() {
		return protein;
	}

	public void setProtein(Protein protein) {
		this.protein = protein;
	}

	public SortedSet<String> getVectorSpaceModel() {
		return vectorSpaceModel;
	}

	public void setVectorSpaceModel(SortedSet<String> vectorSpaceModel) {
		this.vectorSpaceModel = vectorSpaceModel;
	}

	public Map<String, Double> getCumulativeTokenDomainSimilarityScores() {
		return cumulativeTokenDomainSimilarityScores;
	}

	public void setCumulativeTokenDomainSimilarityScores(
			Map<String, Double> cumulativeTokenDomainSimilarityScores) {
		this.cumulativeTokenDomainSimilarityScores = cumulativeTokenDomainSimilarityScores;
	}

	public Double getTotalTokenDomainSimilarityScore() {
		return totalTokenDomainSimilarityScore;
	}

	public void setTotalTokenDomainSimilarityScore(
			Double totalTokenDomainSimilarityScore) {
		this.totalTokenDomainSimilarityScore = totalTokenDomainSimilarityScore;
	}

	public static Map<String, Set<String>> getBlastResultAccessionsToInterproIds() {
		return blastResultAccessionsToInterproIds;
	}

	public static void setBlastResultAccessionsToInterproIds(
			Map<String, Set<String>> blastResultAccessionsToInterproIds) {
		DomainScoreCalculator.blastResultAccessionsToInterproIds = blastResultAccessionsToInterproIds;
	}

	public static Map<String, Set<String>> getBlastResultAccessionsToPfamIds() {
		return blastResultAccessionsToPfamIds;
	}

	public static void setBlastResultAccessionsToPfamIds(
			Map<String, Set<String>> blastResultAccessionsToPfamIds) {
		DomainScoreCalculator.blastResultAccessionsToPfamIds = blastResultAccessionsToPfamIds;
	}
}
