package ahrd.exception;

public class UniprotWebServiceAccessException extends Exception {

	private static final long serialVersionUID = -8658270734576204580L;

	public UniprotWebServiceAccessException() {
		super();
	}

	public UniprotWebServiceAccessException(String m) {
		super(m);
	}

	public UniprotWebServiceAccessException(Exception e) {
		super(e);
	}
}
