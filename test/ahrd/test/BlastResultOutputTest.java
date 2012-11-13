package ahrd.test;

import static org.junit.Assert.assertTrue;
import static ahrd.controller.Settings.getSettings;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.BlastResult;
import ahrd.view.BlastResultOutput;

public class BlastResultOutputTest {

	@Before
	public void setup() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testConstructor() {
		BlastResult br = new BlastResult("accession", 0.0000001,
				"human readable description", 10, 200, 25, 150, 300, 200,
				"swissprot");

		// Most simple output:
		constructBlastResultOutput(br);

		// Output with scores:
		getSettings().setWriteScoresToOutput(true);
		// no description score
		constructBlastResultOutput(br);
		// and with score:
		br.setDescriptionScore(4.75);
		constructBlastResultOutput(br);

		// Output with Token-Set(s):
		getSettings().setWriteTokenSetToOutput(true);
		// no tokens set:
		constructBlastResultOutput(br);
		// and with some tokens:
		br.setTokens(new HashSet<String>(Arrays.asList(new String[] { "human",
				"readable", "description" })));
		constructBlastResultOutput(br);

		// Output with Domain Architecture Similarity scores:
		getSettings().setWriteDomainArchitectureSimilarityScoresToOutput(true);
		// no scores:
		constructBlastResultOutput(br);
		// and with score and domain weight vector:
		br.setDomainSimilarityScore(1.0);
		br.setDomainWeights(Arrays.asList(new Double[] { 0.23, 0.87, 3.245,
				0.0, 0.0 }));
		constructBlastResultOutput(br);
	}

	public void constructBlastResultOutput(BlastResult br) {
		try {
			new BlastResultOutput(br);
		} catch (Exception e) {
			assertTrue("BlastResultOutput constructor threw an error:\n" + e,
					false);
		}
	}
}
