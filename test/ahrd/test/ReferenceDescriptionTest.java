package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.Evaluator;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.Protein;
import ahrd.model.ReferenceDescription;
import ahrd.model.TokenScoreCalculator;
import nu.xom.ParsingException;

public class ReferenceDescriptionTest {

	@Test
	public void testParsingOfReferences() throws IOException {
		TestUtils.initTestSettings();
		String fastaEntry = "AT06g1234 Sheep wool growth factor\nRSSPMSRATVDAAPLLASAAASSGTAPMIEISAAEPKRAPKRVSTTPVTPDRPNSSPPNE\nLIVTVWLFGKMMRSHPTVTRFWPTFRPDW";
		ReferenceDescription rd = ReferenceDescription.constructFromFastaEntry(fastaEntry);
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
	public void testSwissprotBatch1ReferenceTokens()
			throws IOException, MissingAccessionException, MissingProteinException, SAXException, ParsingException {
		Evaluator e = new Evaluator("./test/resources/evaluator_filter_references_test.yml");
		e.initializeProteins();
		e.setupReferenceDescriptions();
		for (Iterator<Map.Entry<String, Protein>> iterator = e.getProteins().entrySet().iterator(); iterator
				.hasNext();) {
			Protein p = iterator.next().getValue();
			ReferenceDescription rd = p.getEvaluationScoreCalculator().getReferenceDescription();
			assertTrue("Reference '" + rd.getAccession() + "' has no tokens after AHRD filtering",
					!rd.getTokens().isEmpty());
		}
	}
}
