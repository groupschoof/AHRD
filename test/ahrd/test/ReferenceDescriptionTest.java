package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Set;

import org.junit.Test;

import ahrd.model.ReferenceDescription;
import ahrd.model.TokenScoreCalculator;

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
		Set<String> generatedTokens = TokenScoreCalculator.tokenize(description, new ArrayList<String>());
		assertEquals(7, generatedTokens.size());
		for (String tkn : tokens) {
			assertTrue("Generated tokens do not contain '" + tkn + "'!", generatedTokens.contains(tkn));
		}
	}
}
