package com.hadoop.maxtemperature;

import java.io.IOException;

import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

public class MaxTemperatureMapper
        extends Mapper<LongWritable, Text, Text, DoubleWritable> {

    @Override
    public void map(LongWritable key, Text value, Context context)
            throws IOException, InterruptedException {
        String line = value.toString();
        if (line.indexOf(' ') == -1) {
            return;
        }
        String year = line.substring(0, line.indexOf(' '));
        System.out.println(line);
        Double airTemperature;
        airTemperature = Double.parseDouble(line.substring(line.indexOf(',') + 1));
        context.write(new Text(year), new DoubleWritable(airTemperature));
    }
}