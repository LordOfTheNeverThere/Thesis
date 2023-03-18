package pt.thesis;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.reflect.Array;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;

import org.apache.spark.sql.Column;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.types.ArrayType;

import scala.collection.immutable.ArraySeq;

import static org.apache.spark.sql.functions.col;
import static org.apache.spark.sql.functions.split;
import static org.apache.spark.sql.functions.lit;

public class App {
    public static String PATH = "C:/Users/migue/Documents/Unsynced University/Thesis/data";
    public static SparkSession spark = SparkSession
            .builder().config("spark.master", "local").getOrCreate();

    public static void main(String[] args) {

        // long counts = df.filter(df.col("Correlation").notEqual(1)).count();
        // System.out.println(counts);
        addGroundTruth();

    }

    private static void addGroundTruth() {

        Dataset<Row> corum = spark.read().option("multiline", "true").json("corumPPI.json");

        Dataset<Row> pairwiseCorr = spark.read().format("csv")
                .option("header", "true")
                .load(PATH + "/datasetsTese/BaseModelPairwise.csv");

        pairwiseCorr = pairwiseCorr.select(split(col("_c0"), ";").as("ProteinArray")).drop("_c0");
        pairwiseCorr = pairwiseCorr.withColumn("Corum", lit(0));
        // pairwiseCorr.filter(pairwiseCorr.col("ProteinArray").apply(1).contains("LRP8")).show();
        // This works

        ArrayList<Integer> isInCorumList = new ArrayList<>();
        Dataset<Row> interestingCorumRow = corum.select("subunits(Gene name)");
        pairwiseCorr.withColumn("Corum", when(interestingCorumRow.filter(corum.col("subunits(Gene name)").contains(pairwiseCorr.col("ProteinArray").getItem(0)))));
        pairwiseCorr.foreach((row) -> {

            String proteinA = row.getList(0).get(0).toString();
            String proteinB = row.getList(0).get(1).toString();
            


            if (!interestingCorumRow.filter(corum.col("subunits(Gene name)").contains(proteinA))
                    .filter(corum.col("subunits(Gene name)").contains(proteinB)).isEmpty()) {
                isInCorumList.add(1);
            } else {
                isInCorumList.add(0);
            }
        });
        Dataset<Row> dummyDf = spark.createDataFrame(isInCorumList, null);
        dummyDf.show();


        // foreach((row) -> {
        // System.out.println(row.split(';'));});

    }
}
