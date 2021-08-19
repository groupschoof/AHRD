package ahrd.view;

import java.util.Date;
import java.util.logging.ConsoleHandler;
import java.util.logging.LogRecord;
import java.util.logging.SimpleFormatter;
import java.util.logging.StreamHandler;
import java.util.logging.Level;

public class DualConsoleHandler extends StreamHandler {

    private final ConsoleHandler stderrHandler = new ConsoleHandler();

    public DualConsoleHandler() {
        super(System.out, new SimpleFormatter() {
            private static final String format = "[%1$tF %1$tT] [%2$-7s] %3$s %n";

            @Override
            public synchronized String format(LogRecord lr) {
                return String.format(format,
                        new Date(lr.getMillis()),
                        lr.getLevel().getLocalizedName(),
                        lr.getMessage()
                );
            }
        });
        stderrHandler.setFormatter(new SimpleFormatter() {
            private static final String format = "[%1$tF %1$tT] [%2$-7s] %3$s %n";

            @Override
            public synchronized String format(LogRecord lr) {
                return String.format(format,
                        new Date(lr.getMillis()),
                        lr.getLevel().getLocalizedName(),
                        lr.getMessage()
                );
            }
        });
    }

    @Override
    public void publish(LogRecord record) {
        if (record.getLevel().intValue() <= Level.INFO.intValue()) {
            super.publish(record);
            super.flush();
        } else {
            stderrHandler.publish(record);
            stderrHandler.flush();
        }
    }
}
