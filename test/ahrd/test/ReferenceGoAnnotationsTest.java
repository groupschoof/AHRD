package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import nu.xom.ParsingException;

import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.AHRD;
import ahrd.controller.Utils;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingInterproResultException;
import ahrd.exception.MissingProteinException;
import ahrd.model.Protein;

public class ReferenceGoAnnotationsTest {

	private AHRD ahrd;

	@Before
	public void setUp() throws IOException {
		ahrd = new AHRD(
				"./test/resources/ahrd_input_seq_sim_table_go_prediction.yml");
	}

	@Test
	public void testHasGeneOntologyAnnotations() {
		// Should have GO-Annotations with default test-Settings:
		assertTrue(getSettings().hasGeneOntologyAnnotations());
		getSettings().setPathToGeneOntologyResults(null);
		assertTrue(!getSettings().hasGeneOntologyAnnotations());
		getSettings().setPathToGeneOntologyResults("/not/existing/path.raw");
		assertTrue(!getSettings().hasGeneOntologyAnnotations());
	}

	@Test
	public void testUniqueShortAccessions() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException,
			ParsingException {
		ahrd.setup(false);
		assertNotNull(ahrd.getUniqueBlastResultShortAccessions());
		// Somehow assertEquals does not work on Collections as expected, hence
		// the following work-around:
		List<String> expectedUniqueShortAccessions = Utils
				.fromFile("./test/resources/all_blast_hits_blast8_tabular_searches.txt");
		for (String sa : expectedUniqueShortAccessions) {
			assertTrue(ahrd.getUniqueBlastResultShortAccessions().contains(sa));
		}
		assertEquals(expectedUniqueShortAccessions.size(), ahrd
				.getUniqueBlastResultShortAccessions().size());
		ahrd.getUniqueBlastResultShortAccessions().removeAll(
				expectedUniqueShortAccessions);
		assertEquals(0, ahrd.getUniqueBlastResultShortAccessions().size());
	}

	@Test
	public void testParseReferenceGoAnnotations() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException,
			ParsingException {
		ahrd.setup(false);
		assertNotNull(ahrd.getReferenceGoAnnotations());
		assertTrue(!ahrd.getReferenceGoAnnotations().isEmpty());
		assertEquals(4, ahrd.getReferenceGoAnnotations().size());
		Set<String> refGos = ahrd.getReferenceGoAnnotations()
				.get("AT1G01040.1");
		assertTrue(refGos.contains("GO:0003824"));
		assertTrue(refGos.contains("GO:0003870"));
	}

	@Test
	public void testAnnotatesGoTerms() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException,
			ParsingException, MissingInterproResultException {
		ahrd.setup(false);
		ahrd.assignHumanReadableDescriptions();
		Protein p = ahrd.getProteins().get("gene:chr01.1056:mRNA:chr01.1056");
		assertNotNull(p.getGoResults());
		assertEquals(2, p.getGoResults().size());
		assertTrue(p.getGoResults().contains("GO:0006355"));
		assertTrue(p.getGoResults().contains("GO:0043401"));
	}
}
