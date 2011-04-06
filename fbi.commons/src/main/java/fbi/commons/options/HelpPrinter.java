package fbi.commons.options;

import fbi.commons.TableFormatter;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.spi.Argument;
import org.cyclopsgroup.jcli.spi.Cli;
import org.cyclopsgroup.jcli.spi.Option;
import org.cyclopsgroup.jcli.spi.ParsingContext;

import java.io.PrintWriter;
import java.util.List;

/**
 * Print help messages for ArgumentProcessors
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class HelpPrinter {
    /**
     * The output writer
     */
    public PrintWriter out;

    /**
     * Create a new help printer that prints to std err
     */
    public HelpPrinter() {
        this(new PrintWriter(System.err));
    }

    /**
     * Create a new Help printer
     *
     * @param out the output
     */
    public HelpPrinter(PrintWriter out) {
        this.out = out;
    }

    /**
     * Print the help message
     *
     * @param processor the processor
     */
    public void print(ArgumentProcessor processor){
        ParsingContext context = processor.createParsingContext();

        Cli cli = context.cli();
        Argument argument = context.argument();
        List<Option> options = context.options();

        // print general command information
        out.println("[USAGE]");
        out.print(" " + cli.getName());
        if(cli.getNote() != null && !cli.getNote().isEmpty()){
            out.print(" " + cli.getNote());
        }else{
            // print options
            out.print(" ");
            for (Option option : options) {
                if(!option.isRequired()){
                    out.print("[");
                }
                out.print("-"+option.getName());
                if(!option.isFlag()){
                    out.print(" <" + option.getDisplayName() + ">");
                }

                if(!option.isRequired()){
                    out.print("]");
                }
                out.print(" ");
            }
        }
        out.println();

        // print description
        if(cli.getDescription() != null && !cli.getDescription().isEmpty()){
            out.println();
            out.println("[DESCRIPTION]");
            out.println(cli.getDescription());
        }

        // print options
        if(options.size() > 0){
            out.println();
            out.println("[OPTIONS]");
            TableFormatter table = new TableFormatter(4);
            for (Option option : options) {
                String[] cols = new String[4];
                cols[0] = "-"+option.getName();
                cols[1] = "--"+option.getLongName();
                if(!option.isFlag()){
                    cols[2] = "<"+option.getDisplayName()+">";
                }else{
                    cols[2] = "";
                }
                cols[3] = option.getDescription();
                table.add(cols);
            }
            out.println(table.toString());
        }

        out.flush();
    }
}
