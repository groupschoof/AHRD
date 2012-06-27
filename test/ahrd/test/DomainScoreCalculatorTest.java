package ahrd.test;

import static junit.framework.Assert.*;

import java.io.IOException;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.DomainScoreCalculator;

public class DomainScoreCalculatorTest {

	@Before
	public void initialiseInterproDb() throws IOException {
		TestUtils.initTestSettings();
	}

	@Test
	public void testInitializeBlastResultAccessionsToInterproIds()
			throws IOException {
		DomainScoreCalculator.initializeBlastResultAccessionsToInterproIds();
		assertTrue("IPR012610 should be assigned to DCL2_ARATH",
				DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.get("DCL2_ARATH").contains("IPR012610"));
		assertTrue("IPR012610 should be assigned to DCL2A_ORYSJ",
				DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.get("DCL2A_ORYSJ").contains("IPR012610"));
		assertTrue("IPR020139 should be assigned to DCL2_ARATH",
				DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
						.get("DCL2_ARATH").contains("IPR020139"));
	}
}
