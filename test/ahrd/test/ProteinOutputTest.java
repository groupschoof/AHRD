package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.BlastResult;
import ahrd.model.EvaluationScoreCalculator;
import ahrd.model.Protein;
import ahrd.model.ReferenceDescription;
import ahrd.view.ProteinOutput;

public class ProteinOutputTest {

	private String path2ReferencesFasta;

	@Before
	public void setup() throws IOException {
		TestUtils.initTestSettings();
		// Most simple AHRD settings
		this.path2ReferencesFasta = getSettings().getPathToReferencesFasta();
		getSettings().setPathToReferencesFasta(null);
		getSettings().setWriteBestBlastHitsToOutput(false);
		getSettings().setWriteDomainArchitectureSimilarityScoresToOutput(false);
		getSettings().setWriteScoresToOutput(false);
		getSettings().setWriteTokenSetToOutput(false);
	}

	@Test
	public void testConstructor() {
		Protein prot = new Protein("P727272", 200);
		List<BlastResult> brs1 = TestUtils.mockBlastResultsForDescCalcTest();
		prot.getBlastResults().put("swissprot", brs1);
		List<BlastResult> brs2 = TestUtils.mockTrEmblBlastResults();
		prot.getBlastResults().put("trembl", brs2);

		// Mock AHRD result:
		BlastResult hsbr = brs1.get(0);
		hsbr.setDescriptionScore(4.5);
		prot.getDescriptionScoreCalculator().setHighestScoringBlastResult(hsbr);
		prot.getDescriptionScoreCalculator().setDescriptionHighScore(4.5);

		// Most simple output:
		constructProteinOutput(prot);
		// ...restore pathToReferencesFasta:
		getSettings().setPathToReferencesFasta(this.path2ReferencesFasta);

		// Output with scores:
		ReferenceDescription rd = new ReferenceDescription();
		rd.setAccession("P727272");
		rd.setDescription("Sweet sheep no goat protein");
		rd.setTokens(new HashSet<String>(Arrays.asList(new String[] { "sweet",
				"sheep", "no", "goat", "protein" })));
		prot.setEvaluationScoreCalculator(new EvaluationScoreCalculator(prot));
		prot.getEvaluationScoreCalculator().setReferenceDescription(rd);
		prot.getEvaluationScoreCalculator().setEvalutionScore(0.85);
		prot.getEvaluationScoreCalculator().setEvalScoreMinBestCompScore(0.3);
		prot.getEvaluationScoreCalculator().setTruePositivesRate(0.85);
		prot.getEvaluationScoreCalculator().setFalsePositivesRate(0.0);
		getSettings().setWriteScoresToOutput(true);
		// no description score
		constructProteinOutput(prot);
		// and with score:
		prot.getTokenScoreCalculator().setTokenHighScore(1.5);

		constructProteinOutput(prot);

		// Output with Token-Set(s):
		getSettings().setWriteTokenSetToOutput(true);
		constructProteinOutput(prot);

		// Output with Domain Architecture Similarity scores:
		getSettings().setWriteDomainArchitectureSimilarityScoresToOutput(true);
		// no scores:
		constructProteinOutput(prot);
		// and with score and domain weight vector:
		hsbr.setDomainSimilarityScore(1.0);
		List<Double> dasVect = Arrays.asList(new Double[] { 0.23, 0.87, 3.245,
				0.0, 0.0 });
		prot.setDomainWeights(dasVect);
		hsbr.setDomainWeights(dasVect);
		constructProteinOutput(prot);
	}

	public void constructProteinOutput(Protein prot) {
		try {
			new ProteinOutput(prot);
		} catch (Exception e) {
			assertTrue("ProteinOutput constructor threw an error:\n" + e, false);
		}
	}
}
