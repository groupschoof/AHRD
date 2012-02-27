package ahrd.test;

import static org.junit.Assert.*;
import org.junit.Test;
import ahrd.model.Protein;
import ahrd.model.GeneOntologyResult;
import java.util.Map;
import java.io.IOException;

public class GeneOntologyResultTest {

	public GeneOntologyResultTest() {
		super();
	}

	@Test
	public void testParseGoResult() throws IOException {
		TestUtils.initTestSettings();
		Map<String, Protein> proteinDb = TestUtils.mockProteinDb();
		try {
			GeneOntologyResult.parseGeneOntologyResult(proteinDb);
		} catch (IOException e) {
			e.printStackTrace(System.out);
		}
		assertEquals(1, proteinDb.get("gene:chr01.502:mRNA:chr01.502")
				.getGoResults().size());
		assertEquals(2, proteinDb.get("gene:chr01.1056:mRNA:chr01.1056")
				.getGoResults().size());
		assertEquals("ribosomal chaperone activity",
				((GeneOntologyResult) proteinDb.get(
						"gene:chr01.502:mRNA:chr01.502").getGoResults()
						.toArray()[0]).getName());
		assertEquals(0.564738, ((GeneOntologyResult) proteinDb.get(
				"gene:chr01.502:mRNA:chr01.502").getGoResults().toArray()[0])
				.getProbability(), 0);
	}

}
