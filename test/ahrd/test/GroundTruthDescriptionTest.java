package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import ahrd.controller.Evaluator;
import ahrd.exception.MissingAccessionException;
import ahrd.model.Protein;
import ahrd.model.GroundTruthDescription;
import ahrd.model.TokenScoreCalculator;

public class GroundTruthDescriptionTest {

	@Test
	public void testParsingOfGroundTruth() throws IOException, MissingAccessionException {
		TestUtils.initTestSettings();
		String fastaEntry = ">AT06g1234 Sheep wool growth factor\nRSSPMSRATVDAAPLLASAAASSGTAPMIEISAAEPKRAPKRVSTTPVTPDRPNSSPPNE\nLIVTVWLFGKMMRSHPTVTRFWPTFRPDW";
		GroundTruthDescription rd = GroundTruthDescription.constructFromFastaEntry(fastaEntry);
		assertEquals("AT06g1234", rd.getAccession());
		assertEquals("Sheep wool growth factor", rd.getDescription());
		assertEquals(4, rd.getTokens().size());
	}

	@Test
	public void testTokenizeDescription() {
		String description = "Sheep wool growth factor putative subfamiLy aCtiviTy";
		String[] tokens = { "sheep", "wool", "growth", "factor", "putative", "subfamily", "activity" };
		Set<String> generatedTokens = TokenScoreCalculator.tokenize(description, new HashSet<String>());
		assertEquals(7, generatedTokens.size());
		for (String tkn : tokens) {
			assertTrue("Generated tokens do not contain '" + tkn + "'!", generatedTokens.contains(tkn));
		}
	}

	@Test
	public void testSwissprotBatch1GroundTruthTokens()
			throws IOException, MissingAccessionException {
		Evaluator e = new Evaluator("./test/resources/evaluator_filter_ground_truth_test.yml");
		e.initializeProteins();
		e.setupGroundTruthDescriptions();
		for (Iterator<Map.Entry<String, Protein>> iterator = e.getProteins().entrySet().iterator(); iterator
				.hasNext();) {
			Protein p = iterator.next().getValue();
			GroundTruthDescription rd = p.getEvaluationScoreCalculator().getGroundTruthDescription();
			assertTrue("Reference '" + rd.getAccession() + "' has no tokens after AHRD filtering",
					!rd.getTokens().isEmpty());
		}
	}
}
