package fbi.commons;

import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;

/**
 * Delegates log messages to our {@link Log} implementations
 *
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */
public class JavaLoggerHandler extends Handler{
    @Override
    public void publish(LogRecord record) {
        Level level = record.getLevel();
        if(level.intValue() == Level.INFO.intValue()){
            Log.info(record.getLoggerName(), record.getMessage());
            return;
        }else if(level.intValue() == Level.WARNING.intValue()){
            Log.warn(record.getLoggerName(), record.getMessage());
        }else if(level.intValue() == Level.SEVERE.intValue()){
            Log.error(record.getMessage(), record.getThrown());
        }else {
            Log.debug(record.getLoggerName(), record.getMessage(), record.getThrown());
        }
    }

    @Override
    public void flush() {
    }

    @Override
    public void close() throws SecurityException {
    }
}
