package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import ahrd.controller.Settings;

public class GOdbSQL {

	/**
	 * Joins all Gene Ontology (GO) terms given in argument goTermAccs with
	 * their parental terms in a set of GO term accessions.
	 * 
	 * @param goTermAccs
	 * @param goDB
	 * @return Set<String>
	 */
	public static Set<String> uniqueGOAccessions(Collection<String> goTermAccs, Map<String, GOterm> goDB) {
		Set<String> g = new HashSet<String>();
		for (String a : goTermAccs)
			if (goDB.containsKey(a))
				g.addAll(goDB.get(a).getParentAccessions());
			else // deprecated or simply non-existing GO Term:
				g.add(a);
		return g;
	}

	/**
	 * Connects to the MySQL Gene Ontology database specified in Settings.
	 * 
	 * @return Connection
	 * @throws SQLException
	 */
	public static Connection connectToGeneOntologyDb() throws SQLException {
		return DriverManager.getConnection(getSettings().getGoDbURL(), getSettings().getGoDbUser(),
				getSettings().getGoDbPassword());
	}

	/**
	 * Finds the Gene Ontology (GO) terms that are parent to the argument
	 * 'goTermAccs' and are NOT obsolete.
	 * 
	 * @param goTermAccs
	 * @param includeSelves
	 * @param relationshipTypeId
	 * @param goCon
	 * @return Map<String, GOterm> Keys are the GO accessions in 'goTermAccs'
	 *         and Values are the instances of GOterm parental to them.
	 * @throws SQLException
	 */
	public static Map<String, GOterm> parentGoTermsForAccessions(Collection<String> goTermAccs, Connection goCon)
			throws SQLException {
		Statement stmt = null;
		String query = parentGoTermsForAccessionsSQLQuery(goTermAccs);
		Map<String, GOterm> res = new HashMap<String, GOterm>();
		Map<String, Set<String>> parentGos = new HashMap<String, Set<String>>();
		try {
			stmt = goCon.createStatement();
			ResultSet rs = null;
			boolean hasRows = stmt.execute(query);
			if (hasRows) {
				rs = stmt.getResultSet();
				while (rs.next()) {
					String gAcc = rs.getString(Settings.GO_DB_TERM_TBL_ACCESSION_KEY);
					if (!res.containsKey(gAcc))
						res.put(gAcc, new GOterm(gAcc, rs.getString(Settings.GO_DB_TERM_TBL_NAME_KEY),
								rs.getString(Settings.GO_DB_TERM_TBL_ONTOLOGY_KEY)));
					String descAcc = rs.getString(Settings.GO_DB_DESCENDANT_TERM_TBL_ACC_KEY);
					if (!parentGos.containsKey(descAcc))
						parentGos.put(descAcc, new HashSet<String>());
					parentGos.get(descAcc).add(gAcc);
				}
			}
		} finally {
			if (stmt != null)
				stmt.close();
		}
		for (String a : parentGos.keySet())
			res.get(a).setParentAccessions(parentGos.get(a));
		return res;
	}

	/**
	 * Concatonates the argument 'goTermAccs' into a single comma separated
	 * String.
	 * 
	 * @param goTermAccs
	 * @return The GO Term accessions concatonated and separated by ","
	 */
	public static String catGoTermAccessions(Collection<String> goTermAccs) {
		StringBuffer gta = new StringBuffer();
		int i = 0;
		for (String g : goTermAccs) {
			gta.append("'" + g + "'");
			if (i < goTermAccs.size() - 1)
				gta.append(",");
			i++;
		}
		return gta.toString();
	}

	/**
	 * Generates a valid SQL query in table term and graph_path to SELECT all
	 * parental Gene Ontology (GO) terms for argument goTermAccs. <i>Note</i>
	 * that the query will also SELECT the GO terms indicated by argument
	 * goTermAccs.
	 * 
	 * @param goTermAccs
	 * @return String
	 */
	public static String parentGoTermsForAccessionsSQLQuery(Collection<String> goTermAccs) {
		return "SELECT t.*, to_root.relation_distance, child.acc as desc_acc "
				+ "FROM graph_path res LEFT JOIN term t ON t.id = res.term1_id "
				+ "LEFT JOIN graph_path to_root ON t.id = to_root.term2_id "
				+ "LEFT JOIN term child ON child.id = res.term2_id " + "WHERE "
				+ "res.term1_id != (SELECT r.id FROM term r WHERE r.is_root = 1) " + "AND child.acc in ( "
				+ goTermAccessionsForAccsOrSynonyms(goTermAccs) + " ) "
				+ "AND to_root.term1_id = (SELECT r.id FROM term r WHERE r.is_root = 1) " + "AND t.is_obsolete = 0 "
				+ "GROUP BY t.id ORDER BY to_root.relation_distance ASC";
	}

	/**
	 * Generates the SQL SELECT Query to obtain the up-to-date GO Term
	 * Accessions even of deprecated GO Term Accessions in 'goTermAccs'.
	 * 
	 * @param goTermAccs
	 * @return String - The valid SQL SELECT statement querying the tables term
	 *         and left joined term_synonym in order to obtain current Go Term
	 *         Accessions
	 */
	public static String goTermAccessionsForAccsOrSynonyms(Collection<String> goTermAccs) {
		String catGoAccs = catGoTermAccessions(goTermAccs);
		return "SElECT t.acc FROM term t LEFT JOIN term_synonym s ON t.id = s.term_id WHERE t.acc in (" + catGoAccs
				+ ") OR s.term_synonym in (" + catGoAccs + ")";
	}
}