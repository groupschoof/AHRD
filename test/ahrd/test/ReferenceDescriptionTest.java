package ahrd.test;

import static ahrd.controller.Settings.setSettings;
import static ahrd.model.ReferenceDescription.constructFromFastaEntry;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.Evaluator;
import ahrd.controller.Settings;
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
		ReferenceDescription rd = constructFromFastaEntry(fastaEntry);
		assertEquals("AT06g1234", rd.getAccession());
		assertEquals("Sheep wool growth factor", rd.getDescription());
		assertEquals(4, rd.getTokens().size());
	}

	@Test
	public void testTokenizeDescription() {
		String description = "Sheep wool growth factor putative subfamiLy aCtiviTy";
		String[] tokens = { "sheep", "wool", "growth", "factor", "putative", "subfamily", "activity" };
		Set<String> generatedTokens = TokenScoreCalculator.tokenize(description, new ArrayList<String>());
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
		e.setupReferences();
		for (Iterator<Map.Entry<String, Protein>> iterator = e.getProteins().entrySet().iterator(); iterator
				.hasNext();) {
			Protein p = iterator.next().getValue();
			ReferenceDescription rd = p.getEvaluationScoreCalculator().getReferenceDescription();
			assertTrue("Reference '" + rd.getAccession() + "' has no tokens after AHRD filtering",
					!rd.getTokens().isEmpty());
		}
	}

	@Test
	public void testReferenceDescriptionBlacklistingAndFiltering() throws IOException {
		setSettings(new Settings(
				Paths.get("test", "resources", "evaluator_filter_references_example_input.yml").toString()));
		ReferenceDescription rd1 = constructFromFastaEntry(
				"ATMG00450.1 hypothetical protein\nMVVTAYPKSSAGMGVTVLPEYLKQSSYEAYSRPYSAFFLSGCTKQERSPLLARRLVDAWL");
		ReferenceDescription rd2 = constructFromFastaEntry(
				"AT1G31870.1 unknown protein\nMAGNQSLKDYLKKYESSDVVEKKKKKKKQKKPSKPEPRGVLVVDEDPVWQKQVDPEEDEN");
		ReferenceDescription rd3 = constructFromFastaEntry(
				"AT1G75110.1 REDUCED RESIDUAL ARABINOSE 2\nMAGRRDRIQQLRGSRIAIAIFVGILIGCVCSVLFPNGFFNSGSSLIANEERISKSTSTDG");
		ReferenceDescription rd4 = constructFromFastaEntry(
				"AT1G75080.1 BRASSINAZOLE-RESISTANT 1; DNA binding / transcription regulator/ transcription repressor\nMTSDGATSTSAAAAAAAAAAARRKPSWRERENNRRRERRRRAVAAKIYTGLRAQGDYNLP");
		assertNotNull("Expected ReferenceDescription rd1 but got null.", rd1);
		assertNotNull("Expected ReferenceDescription rd2 but got null.", rd2);
		assertNotNull("Expected ReferenceDescription rd3 but got null.", rd3);
		assertNotNull("Expected ReferenceDescription rd4 but got null.", rd4);
		assertTrue("Expected rd1 to not have any tokens, instead has: (" + rd1.getTokens() + ")",
				rd1.getTokens().isEmpty());
		assertTrue("Expected rd2 to not have any tokens, instead has: (" + rd2.getTokens() + ")",
				rd2.getTokens().isEmpty());
		assertTrue("Expected ReferenceDescription to contain Token 'reduced', instead has: (" + rd3.getTokens() + ").",
				rd3.getTokens().contains("reduced"));
		assertTrue("Expected ReferenceDescription to contain Token 'residual', instead has: (" + rd3.getTokens() + ").",
				rd3.getTokens().contains("residual"));
		assertTrue(
				"Expected ReferenceDescription to contain Token 'arabinose', instead has: (" + rd3.getTokens() + ").",
				rd3.getTokens().contains("arabinose"));
		assertTrue("Expected ReferenceDescription to contain Token '2', instead has: (" + rd3.getTokens() + ").",
				rd3.getTokens().contains("2"));
		assertTrue("Expected ReferenceDescription to contain Token 'brassinazole', instead has: (" + rd4.getTokens()
				+ ").", rd4.getTokens().contains("brassinazole"));
		assertTrue(
				"Expected ReferenceDescription to contain Token 'resistant', instead has: (" + rd4.getTokens() + ").",
				rd4.getTokens().contains("resistant"));
		assertTrue("Expected ReferenceDescription to contain Token '1', instead has: (" + rd4.getTokens() + ").",
				rd4.getTokens().contains("1"));
		assertTrue("Expected ReferenceDescription to contain Token 'dna', instead has: (" + rd4.getTokens() + ").",
				rd4.getTokens().contains("dna"));
		assertTrue("Expected ReferenceDescription to contain Token 'binding', instead has: (" + rd4.getTokens() + ").",
				rd4.getTokens().contains("binding"));
		assertTrue("Expected ReferenceDescription to contain Token 'transcription', instead has: (" + rd4.getTokens()
				+ ").", rd4.getTokens().contains("transcription"));
		assertTrue(
				"Expected ReferenceDescription to contain Token 'regulator', instead has: (" + rd4.getTokens() + ").",
				rd4.getTokens().contains("regulator"));
		assertTrue("Expected ReferenceDescription to contain Token 'transcription', instead has: (" + rd4.getTokens()
				+ ").", rd4.getTokens().contains("transcription"));
		assertTrue(
				"Expected ReferenceDescription to contain Token 'repressor', instead has: (" + rd4.getTokens() + ").",
				rd4.getTokens().contains("repressor"));
	}
}
