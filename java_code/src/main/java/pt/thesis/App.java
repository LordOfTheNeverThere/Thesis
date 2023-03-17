package pt.thesis;

import tech.tablesaw.api.Table;
import tech.tablesaw.io.csv.*;
public class App 
{
    public static void main( String[] args )
    {
    
        Table t = Table.read().csv(CsvReadOptions.builder("../data/datasetsTese/BaseModelPairwise.csv"));
        System.out.println(t.print());
    }
}
