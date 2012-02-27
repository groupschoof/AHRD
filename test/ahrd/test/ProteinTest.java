package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Map;

import org.junit.Test;

import ahrd.exception.MissingAccessionException;
import ahrd.model.Protein;

public class ProteinTest {

	public ProteinTest() {
		super();
	}

	@Test
	public void testProtConstrInitsTrainScorCalcWhenNeeded() throws IOException {
		TestUtils.initTestSettings();
		Protein p = TestUtils.mockProtein();
		// The above Test-Settings are set to Evaluator-Mode and set to Write
		// out
		// the best Blast-Hits, so we should have a non-null
		// EvaluationScoreCalculator:
		assertTrue(
				"Expected mocked protein's EvaluationScoreCalculator not te be null!",
				p.getEvaluationScoreCalculator() != null);
		// Test the opposite:
		getSettings().setWriteBestBlastHitsToOutput(false);
		getSettings().setPathToReferencesFasta(null);
		p = TestUtils.mockProtein();
		assertTrue(
				"Expected mocked protein's EvaluationScoreCalculator te be null!",
				p.getEvaluationScoreCalculator() == null);
	}

	@Test
	public void testProteinConstructionFromFasta()
			throws MissingAccessionException {
		// NOTE: We expect this method to be called with a split String, read
		// from a FASTA-File.
		// The splitting is done for '>'!
		String input = "MySequence | more information | we don't need\n"
				+ "MADDSKFCFFLVSTFLLLAVVVNVTLAANYVPGDDILLNCGGPDNLPDADGRKWGTDIGS\n"
				+ "KYMLGSKSSTSDAADQKPSVPQVPFMSARVFQSEFTYSFPVASGRKFVRLYFYPSSYNKL\n"
				+ "NATNAIFNVASGQYTLLRNFSAAQTAEALNYDYLTKEFSINVRWKYLHCNVISPSGDTGM\n"
				+ "SPGYDASMTDSRSSGISMSIGGRSLASEDSDGLTPSAVFSQIMNPKGR*";
		String sequence = "MADDSKFCFFLVSTFLLLAVVVNVTLAANYVPGDDILLNCGGPDNLPDADGRKWGTDIGS"
				+ "KYMLGSKSSTSDAADQKPSVPQVPFMSARVFQSEFTYSFPVASGRKFVRLYFYPSSYNKL"
				+ "NATNAIFNVASGQYTLLRNFSAAQTAEALNYDYLTKEFSINVRWKYLHCNVISPSGDTGM"
				+ "SPGYDASMTDSRSSGISMSIGGRSLASEDSDGLTPSAVFSQIMNPKGR*";
		Protein prot = Protein.constructFromFastaEntry(input);
		assertNotNull(
				"Protein should be constructable from a single FASTA-Entry.",
				prot);
		assertEquals(
				"Constructed protein should have it's Accession correctly set",
				"MySequence", prot.getAccession());
		assertEquals(((new Integer(sequence.length())).longValue()),
				(prot.getSequenceLength().longValue()));
	}

	@Test
	public void testInitialisationOfProteinsFromFasta() throws IOException {
		TestUtils.initTestSettings();
		Map<String, Protein> prot_db = null;
		try {
			prot_db = Protein.initializeProteins(getSettings()
					.getProteinsFasta());
		} catch (Exception e) {
			e.printStackTrace(System.out);
		}
		assertTrue(prot_db instanceof Map);
		assertNotNull(prot_db);
		assertEquals(2, prot_db.size());
	}

}
