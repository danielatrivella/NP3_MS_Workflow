#!/usr/bin/env node

/**
 * Module dependencies.
 */

var shell = require('shelljs');
shell.config.silent = true;
var program = require('commander');
//var child_process = require('child_process')

/**
 *  Auxiliary functions
 */

/**
 * @param t process.hrtim(t_initial) - a time diff
 * @returns {string}
 */
function printTimeElapsed(t) {
    return ('Time Elapsed: ' + (t[0]/60).toFixed(0) + "min " +
        (t[0]%60).toFixed(0) + "s " + (t[1]/1000000).toFixed(2) + "ms") ;
}

function printTimeElapsed_bigint(t_start,t_end) {
    const t_elapsed = Number(t_end - t_start)/1000000000; // t is in nanoseconds, compute the diff and convert to seconds
    return 'Time Elapsed: ' + (t_elapsed/60).toFixed(2) + 'min';
}

function parseDecimal(val) {
    return parseInt(val, 10);
}

function toupper(val_txt) {
    return val_txt.toUpperCase();
}

function checkJobNameMaxLength(job_name) {
    if (job_name.length > 80) {
        console.error('\nERROR. The job name is too big, it must have less than 80 characters. ' +
            'Number of characters in the provided job name: '+job_name.length);
        process.exit(1);
    }
    return job_name;
}

function parseBFLAGcutoff(val) {
    val = toupper(val);
    if (val === "FALSE") {
        return val;
    }
    var val_float = parseFloat(val);
    if (isNaN(val_float) || val_float < 0.0) {
        console.error('\nERROR. Wrong bflag_cutoff parameter value. It should be a positive numeric value.');
        process.exit(1);
    }
    return val_float;
}

function parseNOISEcutoff(val) {
    val = toupper(val);
    if (val === "FALSE") {
        return val;
    }
    var val_float = parseFloat(val);
    if (isNaN(val_float) || val_float < 0.0) {
        console.error('\nERROR. Wrong noise_cutoff parameter value. It should be a positive numeric value.');
        process.exit(1);
    }
    return val_float;
}


function splitListFloat(val) {
    var floatList = val.split(",").map(parseFloat);

    if (floatList.some(isNaN) || floatList.length < 2)
    {
        console.error('\nERROR. Wrong rt_tolerance parameter value \'' + val+
            '\'. The retention time tolerance must be two numeric values separated by a comma (e.g. \'1,2\'). Execution aborted.');
        process.exit(1);
    }

    return floatList;
}

function splitList(val) {
    return val.split(",");
}

function increaseVerbosity(v, total) {
    return total + 1;
}

function basename(str) {
    // remove trailing separator
    while (str[str.length - 1] === "/" || str[str.length - 1] === "\\") {
        str = str.substring(0, str.length - 1);
    }
    index_sep = str.lastIndexOf("/");
    if (str.lastIndexOf("\\") > index_sep) {
        index_sep = str.lastIndexOf("\\")
    }

    return str.substr(index_sep + 1);
}


function basedir(str) {
    // remove trailing separator
    while (str[str.length - 1] === "/" || str[str.length - 1] === "\\") {
        str = str.substring(0, str.length - 1);
    }
    index_sep = str.lastIndexOf("/");
    if (str.lastIndexOf("\\") > index_sep) {
        index_sep = str.lastIndexOf("\\")
    }

    return str.substr(0, index_sep+1);
}

function convertIonMode(mode) {
    mode = parseDecimal(mode);

    if (![1,2].includes(mode)) {
        console.error('\nERROR. Wrong ion_mode parameter value. The ion mode must be a positive numeric value in {1,2}. Execution aborted.');
        process.exit(1);
    }

    if (mode === 2)
        return(-1);
    else
        return(mode);
}

function convertMethodCorr(method) {
    if (!["pearson", "kendall", "spearman"].includes(method)) {
        console.error('\nERROR. Wrong method parameter value. The correlation method must be one of {"pearson", "kendall", "spearman"}. Execution aborted.');
        process.exit(1);
    }

    return(method);
}

function convertSimilarityFunc(simFunc) {
    if (!["np3_shifted_cosine", "spec2vec"].includes(simFunc)) {
        console.error('\nERROR. Wrong similarity function parameter value. The similarity function must be one of {"np3_shifted_cosine", "spec2vec"}. Execution aborted.');
        process.exit(1);
    }

    return(simFunc);
}

function callPlotBasePeakIntDistribution(path_clustering_count, bflag_cutoff_factor, logOutputPath, verbose)
{
    const start_plot = process.hrtime.bigint();
    // convert the bflag cutoff to zero when it is disabled
    if (bflag_cutoff_factor === "FALSE") {
        bflag_cutoff_factor = 0;
    }
    var step_name = '*** Plotting the base peak intensity distribution of the clustering counts ***\n';
    console.log(step_name);
    var resExec = shell.exec(python3()+' '+__dirname+'/src/plot_basePeakInt_distribution.py ' + path_clustering_count +
            ' ' + bflag_cutoff_factor, {async: false, silent: (verbose <= 0)});

    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR');
        shell.ShellString(step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(logOutputPath);
    } else {
        var done_msg = '\nDONE! '+printTimeElapsed_bigint(start_plot, process.hrtime.bigint())+"\n";
        console.log(done_msg);
        shell.ShellString(step_name +
            resExec.stdout+'\n'+resExec.stderr + done_msg).toEnd(logOutputPath);
    }
}

function callClustering(options, output_path, specs_path) {
    const start_clust = process.hrtime.bigint();
    console.log('\n*** Step 3 - Clustering the pre-processed MS2 spectra using the MS1 information *** \n');
    // copy MSCluster model
    //shell.cp(options.model_dir+"/LTQ_TRYP_config.txt", output_path);
    // copy metadata table
    shell.cp(options.metadata, output_path);
    processed_dir = options.raw_data_path+'/'+options.processed_data_name;
    var logOutputPath;

    // call not blank samples clustering, SAMPLE_TYPE != blank
    if (shell.test('-e', specs_path+"/data_lists") && !(shell.ls("-A", specs_path+"/data_lists").toString() === ""))
    {
        shell.ls("-A", specs_path+"/data_lists").forEach(function (spec) {

            var out_name = spec.toString().split(".")[0];
            shell.mkdir("-p", output_path+"/outs/"+out_name);
            logOutputPath = callMSCluster(options, options.similarity,specs_path+"/data_lists/"+spec.toString(),
                out_name, output_path, 0, false, 1, 2);
            callCountSpectraBySample(output_path+"/outs/"+out_name, out_name, options.metadata,
                0, processed_dir, options.mz_tolerance, logOutputPath,
                options.verbose);
        });
    }

    // call blank samples clustering, SAMPLE_TYPE == blank
    if (shell.test('-e', specs_path+"/blank_lists") && !(shell.ls("-A", specs_path+"/blank_lists").toString() === ""))
    {
        shell.ls("-A", specs_path+"/blank_lists").forEach(function (spec) {

            var out_name = spec.toString().split(".")[0];
            shell.mkdir("-p", output_path+"/outs/"+out_name);
            logOutputPath = callMSCluster(options, options.similarity_blank, specs_path+"/blank_lists/"+spec.toString(),
                out_name, output_path, 0, false, 1, 2);
            callCountSpectraBySample(output_path+"/outs/"+out_name, out_name, options.metadata,
                0, processed_dir, options.mz_tolerance, logOutputPath,
                options.verbose);
        });
    }

    // call data collection integration step clustering
    if (shell.test('-e', specs_path+"/batch_lists") && !(shell.ls("-A", specs_path+"/batch_lists").toString() === ""))
    {
        shell.ls("-A", specs_path+"/batch_lists").forEach(function (spec) {

            var out_name = spec.toString().split(".")[0];
            shell.mkdir("-p", output_path+"/outs/"+out_name);
            logOutputPath = callMSCluster(options, options.similarity,specs_path+"/batch_lists/"+spec.toString(),
                out_name, output_path, 0, false, 1, 2);

            if (shell.test("-e", output_path+'/outs/'+out_name+'_0') || shell.test("-e", output_path+'/outs/'+out_name+'_1'))
            {
                callCountSpectraBySubClusterID(output_path+"/outs/", out_name, 0, options.metadata,
                    processed_dir, options.mz_tolerance,  logOutputPath,
                    options.verbose);
            } else {
                callCountSpectraBySample(output_path+"/outs/"+out_name, out_name, options.metadata,
                    0, processed_dir, options.mz_tolerance, logOutputPath,
                    options.verbose);
            }
        });
    }

    console.log('*** Integrating batches ***\n');
    shell.mkdir("-p", output_path+"/outs/"+options.output_name);
    logOutputPath = callMSCluster(options, options.similarity,specs_path+'/out_spec_lists.txt', options.output_name, output_path,
        1, true, options.min_peaks_output, 1);

    if (shell.test("-e", output_path+'/outs/B_1') || shell.test("-e", output_path+'/outs/B_1_0') || shell.test("-e", output_path+'/outs/B_1_1'))
    {
        callCountSpectraBySubClusterID(output_path+"/outs/", options.output_name, 1,
            options.metadata, processed_dir, options.mz_tolerance, logOutputPath,
            options.verbose);
    } else {
        callCountSpectraBySample(output_path+"/outs/"+options.output_name, options.output_name,
            options.metadata, 1, processed_dir, options.mz_tolerance,
            logOutputPath, options.verbose);
    }
    var step_time = '\nDONE clustering and counting spectra!\n' + printTimeElapsed_bigint(start_clust, process.hrtime.bigint())+ "\n";
    console.log(step_time);
    shell.ShellString(step_time).toEnd(logOutputPath);

    // for analysing the clustering counts
    callAnalyseCount(output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_spectra.csv",
        output_path+"/outs/"+options.output_name+"/count_tables/analyseCountClustering",
        logOutputPath);

    // concatenate spectra peak list
    callExtractPeakList(options.output_name, output_path+"/outs/"+options.output_name+"/mgf/",
        output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_peak_area.csv",
        output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_spectra.csv",
        options.fragment_tolerance, options.scale_factor, logOutputPath);

    // Create Groups coluns listed at metadata
    callGroupsfunc(options.metadata, 
        output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_peak_area.csv",
        logOutputPath, options.verbose);
    callGroupsfunc(options.metadata, 
        output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_spectra.csv",
        logOutputPath, options.verbose);

    // call plot basePeakInt distribution!!
    callPlotBasePeakIntDistribution(output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_spectra.csv",
        options.bflag_cutoff, logOutputPath, options.verbose);

    // call aggregation of not fragmented MS1 peaks
    // if exists options.raw_data_path+'/'+options.processed_data_name+'/'+"MS1_list_no_MS2.csv"
    if (shell.test('-e', options.raw_data_path+'/'+options.processed_data_name+"/MS1_list_no_MS2.csv"))
        callCleanNoMs2Counts(options.raw_data_path+'/'+options.processed_data_name+"/MS1_list_no_MS2.csv" ,
            options.metadata, output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_peak_area_MS1.csv",
            options.mz_tolerance, options.rt_tolerance[1], options.method, logOutputPath,
            options.verbose);

    step_time = '\n======\nFinish Steps 3 and 4! ' + printTimeElapsed_bigint(start_clust, process.hrtime.bigint())+ "\n======\n";
    console.log(step_time);
    shell.ShellString(step_time).toEnd(logOutputPath);
}


function callMSCluster(parms, sim_tol, spec, name, out_path, rt_tol_i, keep_split_mgf, min_numPeak_output,
                       min_verbose) {
    //shell.cd('NP3_MSCluster');
    console.log('*** Calling NP3_MSCluster for: '+name+' *** \n');
    var resExec;

    if (isWindows())
    {
        resExec = shell.exec('"'+__dirname+'\\NP3_MSCluster\\NP3_MSCluster_bin.exe" --list "'+spec+'" --output-name "'+name+'" ' +
            '--out-dir "'+out_path+'\\outs\\'+name+'" --rt-tolerance '+parms.rt_tolerance[rt_tol_i]+
            ' --fragment-tolerance '+parms.fragment_tolerance+' --window '+parms.mz_tolerance+' --similarity '+sim_tol+
            ' --model-dir "'+parms.model_dir+'" --sqs 0.0 --num-rounds '+parms.num_rounds+' --mixture-prob '+parms.mixture_prob+
            ' --tmp-dir "'+__dirname+'\\NP3_MSCluster\\tmp_'+name+'_rmv"'+' --min-peaks-output '+min_numPeak_output+
            ' --scale-factor '+parms.scale_factor +' --verbose-level 10 --output-mgf --assign-charges --major-increment 100 ' +
            '--output-file-size '+parms.max_chunk_spectra, {async:false, silent:(parms.verbose < min_verbose)});
    } else {
        const tmpDir = "/tmp/NP3_MSCluster/"+Math.random().toString(36).substring(2, 10);
        const msclusterTmpDir = tmpDir+"/tmp_"+name+"_rmv";
        shell.mkdir("-p", msclusterTmpDir);
        resExec = shell.exec(__dirname+'/NP3_MSCluster/NP3_MSCluster_bin --list '+spec+' --output-name '+name+' ' +
            '--out-dir '+out_path+'/outs/'+name+' --rt-tolerance '+parms.rt_tolerance[rt_tol_i]+
            ' --fragment-tolerance '+parms.fragment_tolerance+' --window '+parms.mz_tolerance+' --similarity '+sim_tol+
            ' --model-dir '+parms.model_dir+' --sqs 0.0 --num-rounds '+parms.num_rounds+' --mixture-prob '+parms.mixture_prob+
            ' --tmp-dir '+msclusterTmpDir+' --min-peaks-output '+min_numPeak_output+
            ' --scale-factor '+parms.scale_factor +' --verbose-level 10 --output-mgf --assign-charges --major-increment 100 ' +
            '--output-file-size '+parms.max_chunk_spectra, {async:false, silent:(parms.verbose < min_verbose)});
        shell.rm('-rf', tmpDir);
    }
    // save mscluster log file
    var log_output_path = out_path+'/outs/'+name+'/logClusteringOutput';
    shell.ShellString('\n*** Step 3 - Clustering the pre-processed MS2 spectra using the MS1 information *** \n'+
        '\n*** Calling NP3_MSCluster for: '+name+' *** \n\n'+resExec.stdout+'\n'+
        resExec.stderr).to(log_output_path);

    // remove  file and leave mscluster dir
    shell.rm('-rf', __dirname+'/NP3_MSCluster/tmp_'+name+'_rmv');
    shell.rm('-rf', __dirname+'/NP3_MSCluster/out*');
    //shell.cd('..');

    if (resExec.code) {
        //shell.cd('..');
        if (parms.verbose < min_verbose) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('ERROR');
        shell.ShellString('\nERROR\n').toEnd(log_output_path);
        process.exit(1);
    } else {
        console.log('DONE!\n');
        shell.ShellString('DONE!\n').toEnd(log_output_path);
        // merge mgfs
        shell.cat(shell.cat(out_path+'/outs/'+name+'/'+name+'_0_0_mgf_list.txt').split("\n").filter(String)).to(
            out_path+'/outs/'+name+'/mgf/'+name+'_all.mgf');

        // rm unmerged mgfs when they wont be needed (except last step)
        if (!keep_split_mgf) {
            shell.rm('-rf', out_path + '/outs/' + name + '/mgf/*_[0-9].mgf');
        }
        shell.rm('-rf', out_path+'/outs/'+name+'/'+name+'_0_0_clust_list.txt');
        shell.rm('-rf', out_path+'/outs/'+name+'/'+name+'_0_0_mgf_list.txt');
    }

    return log_output_path;
}

function callCountSpectraBySample(out_path, name, metadata, isFinal, processed_dir, mz_tol, logOutputPath, verbose)
{
    var step_name = '*** Step 4 - Counting peak area and spectra by sample '+name+' *** \n';
    console.log(step_name);
    var resExec = shell.exec('Rscript '+__dirname+'/src/count_clust_samples.R '+out_path+' '+name+' '+metadata+' '+
        isFinal+' '+processed_dir+' '+mz_tol,
        {async:false, silent:(verbose <= 0)});

    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR');
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(logOutputPath);
        process.exit(1);
    } else {
        console.log('\nDONE!\n');
        shell.ShellString('\n'+ step_name +
            resExec.stdout+'\n'+resExec.stderr+'\nDONE!\n').toEnd(logOutputPath);
    }
}

function callCountSpectraBySubClusterID(out_path, name, isFinal, metadata, processed_dir, mz_tol, logOutputPath, verbose)
{
    var step_name = '*** Step 4 - Counting peak area and spectra by batch subclusters '+name+' *** \n';
    console.log(step_name);
    var resExec = shell.exec('Rscript '+__dirname+'/src/count_clust_batches.R '+out_path+' '+name+' '+isFinal+' '+
        metadata+' '+processed_dir+' '+mz_tol, {async:false, silent:(verbose <= 0)});

    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR');
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(logOutputPath);
        process.exit(1);
    } else {
        console.log('\nDONE!\n');
        shell.ShellString('\n'+ step_name +
            resExec.stdout+'\n'+resExec.stderr + '\nDONE!\n').toEnd(logOutputPath);
    }
}

function callComputeCorrelation(metadata, counts, method, bio_cutoff, logOutputPath, verbose)
{
    //# params:
    //   #$1 - Path to the CSV batch metadata file containing filenames, sample codes, data collection batches and blanks;
    //   #$2 - Path to the CSV spectra count file;
    //   #$3 - Path to the batches output folder;
    //   #$4 - output name name;
    //   #$5 - Correlation method, one of: 'pearson' (default), 'kendall', or 'spearman'.
    //   # return NA when div/0, all samples = 0 or return 1.1 when sd equals zero
    var step_name = '*** Step 9 - Computing the correlation between the counts table and the samples bioactivity from '+
        basename(counts)+' *** \n';
    console.log(step_name);
    const start_corr = process.hrtime.bigint();
    var resExec = shell.exec('Rscript '+__dirname+'/src/bioactivity_correlation.R '+metadata+' '+counts+' '+method+' '+
        bio_cutoff,
        {async:false, silent:(verbose <= 0)});

    if (resExec.code) {
        // in case of error show all the emmited msgs
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('ERROR\n');
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(logOutputPath);
    } else {
        var done_msg = '\nDONE! '+printTimeElapsed_bigint(start_corr, process.hrtime.bigint())+"\n";
        console.log(done_msg);
        shell.ShellString('\n'+ step_name + resExec.stdout+'\n'+resExec.stderr+done_msg).toEnd(logOutputPath);
    }
}

function callAnalyseCount(counts, out_path, logOutputPath)
{
    var step_name = '*** Analysing the count of spectra  *** \n';
    console.log(step_name);
    var resExec = shell.exec('Rscript '+__dirname+'/src/analyse_count.R '+counts, {async:false, silent:false});

    if (resExec.code) {
        //console.log(resExec.stdout);
        //console.log(resExec.stderr);
        console.log('\nERROR');
        shell.ShellString(step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(logOutputPath);
    } else {
        console.log('\nDONE!\n');

        // save analyse count file
        shell.ShellString(resExec.stdout).to(out_path);
        shell.ShellString('\n'+step_name+resExec.stdout+'\n'+resExec.stderr+
            '\nDONE!\n').toEnd(logOutputPath);
    }
}

function callExtractPeakList(job_name, mgf_dir, counts_area, counts_spectra, bin_size, scale_factor, logOutputPath)
{
    var step_name = '*** Extracting the fragmented peaks list from the clustered MGF and concatenating it to the counts table *** \n';
    console.log(step_name);
    var resExec = shell.exec('Rscript '+__dirname+'/src/extract_peak_list.R '+job_name+' '+mgf_dir+' '+counts_area+' '+bin_size+' '+
        scale_factor+' '+counts_spectra, {async:false});

    if (resExec.code) {
        console.log(resExec.stdout);
        console.log(resExec.stderr);
        console.log('\nERROR');
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(logOutputPath);
    } else {
        console.log('\nDONE!\n');
        shell.ShellString('\n'+ step_name + resExec.stdout+'\n'+resExec.stderr+'\nDONE!\n').toEnd(logOutputPath);
    }
}

// output_path file name with path to save the result
function callCleanNoMs2Counts(quantification_table_path, metadata_path, output_path, mz_tol, rt_tol, method,
                              logOutputPath, verbose)
{
    const start_cleanMS1 = process.hrtime.bigint();
    var step_name = '*** Cleaning the counts table of not fragmented MS1 peaks  *** \n';
    console.log(step_name);
    var resExec = shell.exec('Rscript '+__dirname+'/src/clean_no_MS2_quantification.R  '+metadata_path+' '+output_path+' '+
        quantification_table_path+' '+mz_tol+' '+rt_tol, {async:false, silent:(verbose === 0)});

    if (resExec.code) {
        console.log(resExec.stdout);
        //console.log(resExec.stderr);
        console.log('\nERROR');
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(logOutputPath);
    } else {
        var done_msg = '\nDONE! '+printTimeElapsed_bigint(start_cleanMS1, process.hrtime.bigint())+"\n";
        console.log(done_msg);
        shell.ShellString('\n' + step_name +
            resExec.stdout+'\n'+resExec.stderr + done_msg).toEnd(logOutputPath);
        // call correlation
        callComputeCorrelation(metadata_path, output_path,
            method, 0, logOutputPath, verbose);
    }
}

function callMergeCounts(output_path, output_name, processed_dir, metadata_path, merge_protonated, method, verbose)
{
    const start_merge = process.hrtime.bigint();
    var step_name = '*** Step 8 - Merging the counts table based on the provided chemical annotations *** \n';
    console.log(step_name);
    annotations_merge = "isotopes,adducts,dimers,multicharges,fragments";
    var merge_output_path = output_path + "/count_tables/merge/";

    var resExec = shell.exec('Rscript '+__dirname+'/src/merge_annotation_counts.R '+output_path+' '+annotations_merge+' '+
        metadata_path+' '+processed_dir+' '+merge_protonated, {async:false, silent:(verbose === 0)});

    if (resExec.code) {
        if (verbose === 0) {
            console.log(resExec.stdout);
            // console.log(resExec.stderr);
        }
        shell.ShellString(step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(merge_output_path+'logMergeOutput');
        console.log('\nERROR');
    } else {
        var done_msg = '\nDONE! '+printTimeElapsed_bigint(start_merge, process.hrtime.bigint())+"\n";
        console.log(done_msg);
        shell.ShellString(step_name +
            resExec.stdout+'\n'+resExec.stderr + done_msg).toEnd(merge_output_path+'logMergeOutput');

        // run corr if the metadata was provided and at least one symbolic cluster was created
        if (typeof metadata_path != "undefined" && shell.test('-e',merge_output_path +
            output_name + '_peak_area_merged_ann.csv')) {
            // call groups
            callGroupsfunc(metadata_path,
                merge_output_path + output_name + '_peak_area_merged_ann.csv',
                merge_output_path+'logMergeOutput', verbose);

            // call correlation
            callComputeCorrelation(metadata_path, merge_output_path +
                output_name + '_peak_area_merged_ann.csv',
                method, 0, merge_output_path+'logMergeOutput', verbose);
            // call correlation
            callComputeCorrelation(metadata_path, merge_output_path +
                output_name + '_spectra_merged_ann.csv',
                method, 0, merge_output_path+'logMergeOutput', verbose);
        }
        var done_msg = '\n======\nFinish Step 8 '+printTimeElapsed_bigint(start_merge, process.hrtime.bigint())+"\n======\n";
        console.log(done_msg);
        shell.ShellString(done_msg).toEnd(merge_output_path+'logMergeOutput');
    }
}

// call the clean step
function callCleanClusteringCounts(parms, output_path, mz_tol, rt_tol, bin_size, processed_dir,
                                   log_clustered_spec_comparison)
{
    // console.log('outp '+ output_path+ " rt "+rt_tol+" bin "+ bin_size)
    const start_clean = process.hrtime.bigint();
    var step_name = '*** Step 5 - Cleaning the clustering counts tables *** \n';
    console.log(step_name);
    if (processed_dir === "") {
        processed_dir = parms.raw_data_path + '/' + parms.processed_data_name;
    }

    var clean_output_path = output_path+"/count_tables/clean/";

    var resExec = shell.exec('Rscript '+__dirname+'/src/clean_spectra_quantification.R '+parms.metadata+' '+output_path+' '+
        processed_dir+' '+mz_tol+' '+ parms.similarity+' '+ rt_tol+' '+bin_size+' '+
        parms.scale_factor+' '+parms.ion_mode+' '+parms.bflag_cutoff+' '+parms.noise_cutoff+' '+parms.max_chunk_spectra,
        {async:false,
        silent:(parms.verbose === 0)});
    // print the clustered spectra comparision in the output log after the folder creation (counts_table/clean)
    shell.ShellString(log_clustered_spec_comparison + "\n" + step_name).toEnd(clean_output_path+"logCleanOutput");

    if (resExec.code) {
        if (parms.verbose === 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        shell.ShellString(resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(clean_output_path+"logCleanOutput");
        console.log('\nERROR');
    } else {
        var out_done_log = '\nDONE! '+printTimeElapsed_bigint(start_clean, process.hrtime.bigint())+"\n";
        console.log(out_done_log);
        shell.ShellString(resExec.stdout+'\n'+resExec.stderr+out_done_log).toEnd(clean_output_path+"logCleanOutput");

        output_name = basename(output_path);
        // analyse the clean clustering count
        callAnalyseCount(clean_output_path+output_name+"_spectra_clean.csv",
            clean_output_path+"analyseCountClusteringClean",
            clean_output_path+"logCleanOutput");
        // call the pairwise comparision of the clean mgf
        var out_pairComp = callPairwiseComparision(output_name+'_clean',
            output_path+'/molecular_networking/similarity_tables/',
            output_path+'/mgf/'+output_name+'_clean.mgf',
             bin_size, parms.scale_factor, parms.trim_mz, parms.max_shift,
            parms.parallel_cores, parms.similarity_function, parms.verbose);
        shell.ShellString(out_pairComp).toEnd(clean_output_path+"logCleanOutput");

        callGroupsfunc(parms.metadata,
            clean_output_path+output_name+"_peak_area_clean.csv",
            clean_output_path+"logCleanOutput", parms.verbose);
        // TODO Future alignment function will probably be here

        var out_done_log = '\n======\nFinish Step 5 '+printTimeElapsed_bigint(start_clean, process.hrtime.bigint())+' (without the clustered spectra comparison time) \n======\n';
        console.log(out_done_log);
        shell.ShellString(out_done_log).toEnd(clean_output_path+"logCleanOutput");
    }

    return(resExec.code);
}

function callAnnotateCleanCounts(parms, output_path, mz_tol, fragment_tol, rt_tol,
                                 absolute_ms2_int_cutoff)
{
    // console.log('outp '+ output_path+ " rt "+rt_tol+" bin "+ bin_size)
    var step_name = '*** Step 7 - Annotating ionization variants in the clean counts table and creating the molecular network of annotations (IVAMN) *** \n';
    console.log(step_name);
    const start_ann = process.hrtime.bigint();

    var resExec = shell.exec('Rscript '+__dirname+'/src/run_annotate_spectra_molecular_network.R '+parms.metadata+' '+output_path+' '+
        parms.rules+' '+mz_tol+' '+ fragment_tol+' '+ rt_tol+' '+absolute_ms2_int_cutoff+' '+parms.ion_mode+' '+
        parms.scale_factor+' '+parms.max_chunk_spectra, {async:false, silent:(parms.verbose === 0)});

    if (resExec.code) {
        if (parms.verbose === 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        shell.ShellString(step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(output_path+"/count_tables/clean/logAnnotateOutput");
        console.log('\nERROR');
        return(resExec.code);
    } else {
        console.log('\nDONE!\n');
        shell.ShellString(step_name +
            resExec.stdout + '\n' + resExec.stderr + '\nDONE!\n').toEnd(output_path+"/count_tables/clean/logAnnotateOutput");
    }

    var step_name = '*** Step 7 - Assigning the putative [M+H]+ spectra representatives in the IVAMN *** \n';
    console.log(step_name);
    const start_protonated = process.hrtime.bigint();
    var output_name = basename(output_path);
    var mn_annotations_path = output_path + '/molecular_networking/' +
        output_name + "_ivamn.selfloop";
    var peak_area_clean_path = output_path + "/count_tables/clean/" +
        output_name+"_peak_area_clean_ann.csv";

    var resExec = shell.exec(python3()+' '+__dirname+'/src/mn_annotations_assign_protonated_representative.py '+
        mn_annotations_path+' '+peak_area_clean_path, {async:false, silent:(parms.verbose === 0)});

    if (resExec.code) {
        if (parms.verbose === 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR');
        shell.ShellString(step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(output_path+"/count_tables/clean/logAnnotateOutput");
    } else {
        console.log('\nDONE!\n');
        shell.ShellString(step_name +
            resExec.stdout + '\n' + resExec.stderr + '\nDONE! '+printTimeElapsed_bigint(start_protonated, process.hrtime.bigint())+'\n').toEnd(output_path+"/count_tables/clean/logAnnotateOutput");
    }

    var finish_step = '\n======\nFinish Step 7 '+printTimeElapsed_bigint(start_ann, process.hrtime.bigint())+"\n======\n";
    console.log(finish_step);
    shell.ShellString(finish_step).toEnd(output_path+"/count_tables/clean/logAnnotateOutput");
    return(resExec.code);
}

function callPairwiseComparision(out_name, out_path, mgf_path, bin_size, scaling_method, trim_mz, max_shift,
                                 cores_parallel, similarity_function, verbose)
{
    const step_name = '*** Step 5 - Computing the pairwise similarity comparisons of the resulting consensus spectra *** \n';
    console.log(step_name);

    const start_comp = process.hrtime.bigint();
    // select the similarity function to be used, one of NP3 shifted cosine or spec2vec - if spec2vec convert mgf dir to _all.mgf file
    if (similarity_function == "spec2vec")
    {
        console.log("* Using the spec2vec function to compute the spectra comparisons");
        console.log("* The scaling_method is automatically set to the power of 0.5, similar to the used model");
        // check if mgf_path is a directory, if so convert it to the single file named <job_name>_all.mgf
        if (shell.test('-d', mgf_path)) {
            mgf_path = mgf_path +'/'+out_name+'_all.mgf'
        }
        var resExec = shell.exec('python '+__dirname+'/src/pairwise_similarity_spec2vec.py '+out_name+' '+mgf_path+' '+out_path+' '+
            bin_size+' '+trim_mz, {async:false, silent:(verbose <= 0)});
    } else {
        console.log("* Using the NP3 shifted cosine function to compute the spectra comparisons");
        var resExec = shell.exec('Rscript '+__dirname+'/src/pairwise_similarity.R '+out_name+' '+mgf_path+' '+out_path+' '+
            bin_size+' '+scaling_method+' '+trim_mz+' '+max_shift+' '+cores_parallel, {async:false, silent:(verbose <= 0)});
    }

    var output_msg = '';
    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR');
        shell.ShellString(step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(out_path+'/logPairwiseComparisonOutput');
        process.exit(1);
    } else {
        const done_msg = '\nDONE! \n======\nFinish pairwise comparison '+printTimeElapsed_bigint(start_comp, process.hrtime.bigint())+'\n======\n';
        console.log(done_msg);
        output_msg = step_name +
            resExec.stdout + '\n' + resExec.stderr + done_msg;
        shell.ShellString(output_msg).toEnd(out_path+'/logPairwiseComparisonOutput');
    }

    return(output_msg)
}

function callCreatMN(out_path, sim_mn, net_top_k, max_component_size, min_matched_peaks, max_chunk_spectra,
                     blank_expansion, verbose)
{
    var step_name ='*** Step 10 - Creating the Spectra Similarity Molecular Network (SSMN) *** \n';
    console.log(step_name);
    const start_mn = process.hrtime.bigint();
    var resExec = shell.exec('Rscript '+__dirname+'/src/molecular_networking.R '+out_path+' '+sim_mn+' '+
        max_chunk_spectra,
        {async:false, silent:(verbose <= 0)});

    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        shell.ShellString(step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(out_path+'/molecular_networking/logMnOutput');
        console.log('\nERROR');
        return(resExec.code)
    } else {
        console.log('\nDONE!\n');
        shell.ShellString(step_name +
            resExec.stdout + '\n' + resExec.stderr + '\nDONE!\n').toEnd(out_path+'/molecular_networking/logMnOutput');
    }
    // SSMN filtered
    var step_name ='*** Filtering the SSMN - minimum matched peaks, top k neighbours and max component size  *** \n';
    console.log(step_name);
    var mn_file = out_path+'/molecular_networking/'+basename(out_path)+"_ssmn_w_"+
        sim_mn.toString().replace(".", "")+".selfloop";

    resExec = shell.exec(python3()+' '+__dirname+'/src/molecular_network_filtering_library.py '+mn_file+' '+net_top_k+' '+
        max_component_size+ ' '+min_matched_peaks,
        {async:false, silent:(verbose <= 0)});

    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(out_path+'/molecular_networking/logMnOutput');
        console.log('\nERROR');
    } else {
        console.log('\nDONE!\n');
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\n' + resExec.stderr + '\nDONE!\n').toEnd(out_path+'/molecular_networking/logMnOutput');
    }
    // protonated mn - SSMN [M+H]+ and IVAMN [M+H]+
    var step_name ='*** Creating the [M+H]+ networks without blanks - protonated SSMN and IVAMN *** \n';
    console.log(step_name);
    var ivamn_file = out_path+'/molecular_networking/'+basename(out_path)+
        "_ivamn.selfloop";
    var clean_table_file = out_path+'/count_tables/clean/'+basename(out_path)+
        "_peak_area_clean_ann.csv";

    resExec = shell.exec(python3()+' '+__dirname+'/src/create_protonated_ssmn_ivamn.py '+ivamn_file+' '+mn_file+' '+
        clean_table_file+' '+blank_expansion+' '+net_top_k+' '+max_component_size+ ' '+min_matched_peaks,
        {async:false, silent:(verbose <= 0)});

    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(out_path+'/molecular_networking/logMnOutput');
        console.log('\nERROR');
    } else {
        var mn_time = '\nDONE! \n======\nFinish Step 10 '+printTimeElapsed_bigint(start_mn, process.hrtime.bigint())+"\n======\n";
        console.log(mn_time);
        shell.ShellString('\n' + step_name +
            resExec.stdout + '\n' + resExec.stderr + '\nDONE!\n' + mn_time).toEnd(out_path+'/molecular_networking/logMnOutput');
    }

    return(resExec.code)
}

function callPreProcessSuggestion(metadata_path, processed_data_dir, out_path, verbose)
{
    console.log('*** Step 2 Suggestion - Suggesting values for some of the pre process parameters *** \n');
    const start_ppsug = process.hrtime.bigint();
    n_bins = 150;
    var resExec = shell.exec(python3()+' '+__dirname+'/src/pp_pw_hist_suggestions.py '+metadata_path+' '+
        processed_data_dir+'/MS1_list_with_MS2.csv '+processed_data_dir+'/MS1_list_with_MS2_noBlank_peak_width_hist.png '+n_bins,
        {async:false, silent: (verbose <= 0)});

    if (resExec.code) {
        console.log(resExec.stdout);
        console.log(resExec.stderr);
        console.log('\nERROR');
    } else {
        var ppsuggest_time = '\nDONE! '+printTimeElapsed_bigint(start_ppsug, process.hrtime.bigint())+"\n";
        console.log(ppsuggest_time);

        // save the suggestion output to the statistics log and to the output log
        shell.ShellString(resExec.stdout).toEnd(out_path);
        shell.ShellString(resExec.stdout+ppsuggest_time).toEnd(processed_data_dir+'logPreProcessOutput');
        return resExec.stdout
    }
    return ''
}

function callGroupsfunc(metadata_path, count_tables, logOutputPath, verbose)
{
    var step_name = '*** Creating groups based on the metadata grouping information *** \n';
    console.log(step_name);
    var resExec = shell.exec(python3()+' '+__dirname+'/src/groups.py --metadata '+metadata_path+' --count_file_path '+count_tables+
    ' -q True', {async:false, silent: (verbose <= 0)});
    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR');
        shell.ShellString('\n' + step_name+
            resExec.stdout + '\nSTDERR:\n' + resExec.stderr + '\nERROR').toEnd(logOutputPath);
    } else {
        console.log('\nDONE!\n');
        shell.ShellString('\n' + step_name +
            resExec.stdout+'\nDONE!\n').toEnd(logOutputPath);
    }

    return ''
}

function callPreProcessData(job, metadata, raw_dir, parms, verbose)
{
    var resExec;
    console.log('*** Step 2 - Pre-processing the raw LC-MS/MS data and enriching the MS2 data with MS1 peak information *** \n');
    const start_pp = process.hrtime.bigint();

    if (parms.fragment_tolerance) // called from the run or clustering cmd, auto process
    {
        resExec = shell.exec('Rscript '+__dirname+'/src/tandem_peak_info_align.R ' + job + ' ' + metadata + ' ' + raw_dir+ ' ' +parms.rt_tolerance[0]+ ' ' +
            parms.fragment_tolerance + ' ' + parms.ppm_tolerance + ' ' + parms.ion_mode+' '+ parms.processed_data_name+' '+
            parms.processed_data_overwrite, {async: false, silent: (verbose <= 0), maxBuffer: 1024*1024*1024});
    } else { // called from the process cmd
        resExec = shell.exec('Rscript '+__dirname+'/src/tandem_peak_info_align.R ' + job + ' ' + metadata + ' ' + raw_dir + ' ' +parms.rt_tolerance+ ' ' +
            parms.mz_tolerance+ ' ' + parms.ppm_tolerance + ' ' + parms.ion_mode+' '+ parms.processed_data_name +' '+
            parms.processed_data_overwrite+ ' ' +parms.peak_width+ ' ' +parms.snthresh+ ' ' +parms.pre_filter+ ' ' +
            parms.min_fraction+ ' ' +parms.bw+ ' '  +parms.bin_size + ' ' +parms.max_features+ ' ' +parms.noise+ ' ' +parms.mz_center_fun+
            ' '+ parms.max_samples_batch_align+' '+ parms.integrate_method+' '+ parms.fit_gauss, {async: false, silent: (verbose <= 0), maxBuffer: 1024*1024*1024});
    }

    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log("\nError msg logged to: " + raw_dir+'/logPreProcessError');
        console.log('\nERROR');
        // save log error
        shell.ShellString("ERROR\n\nOUTPUT:"+resExec.stdout+"\n\nERROR:"+resExec.stderr).to(raw_dir+'/logPreProcessError');
        process.exit(1);
    }

    var pre_process_dir = raw_dir+'/'+parms.processed_data_name+'/';

    // save the pre_process output if this is a new result
    var preprocess_new_result = resExec.stdout.match(/All raw LC\-MS\/MS files were already pre processed!/gm);
    if (!preprocess_new_result) {
        // no previous result was used, then save the complete output
        shell.ShellString(resExec.stdout).to(pre_process_dir+'logPreProcessOutput');
    }

    var reg_stats = resExec.stdout.match(/\* Pre-processing - Statistics [\s\S]*/gm);
    // console.log(reg_stats);
    if (!reg_stats) {
        console.log("ERROR in the pre-processing, no statistics was printed! :(");
        shell.ShellString("ERROR\n\nOUTPUT:"+resExec.stdout+"\n\nERROR:"+resExec.stderr).to(raw_dir+'/logPreProcessError');
        shell.ShellString("ERROR\n\nOUTPUT:"+resExec.stdout+"\n\nERROR:"+resExec.stderr).to(pre_process_dir+'logPreProcessOutput');

        process.exit(1);
    }
    reg_stats[0] = '\n' + reg_stats[0].trim() + '\n';
    if (verbose <= 0) {
        console.log(reg_stats[0]);
    }
    var output_stats;
    // check for a warning
    var reg_warning = reg_stats[0].match(/\*\*\*\*\*\*\*\*\*\*\*\n\* WARNING \*/gm);
    if (reg_warning) {
        // save log statistics and warning
        shell.ShellString(reg_stats[0]).to(pre_process_dir+'/logPreProcessStatisticsWarning');
        console.log('DONE! Check the warning for recommended actions to improve the pre processing result.\n');
        // call the suggestion if not done yet
        if (!reg_stats[0].match(/\* Computing the peak width of the MS1 list with MS2 and removing blanks \*/gm)) {
            reg_stats[0] += callPreProcessSuggestion(metadata, pre_process_dir,
                pre_process_dir + '/logPreProcessStatisticsWarning', verbose);
        }
        output_stats = reg_stats[0];
    } else {
        shell.ShellString(reg_stats[0]).to(pre_process_dir+'/logPreProcessStatistics');
        console.log('DONE!\n');
        // call the suggestion if not done yet
        if (!reg_stats[0].match(/\* Computing the peak width of the MS1 list with MS2 and removing blanks \*/gm)) {
            callPreProcessSuggestion(metadata, pre_process_dir,
                pre_process_dir + '/logPreProcessStatistics', verbose);
        }
        output_stats = undefined;
    }
    var pp_end = "\n======\nFinish Step 2 "+printTimeElapsed_bigint(start_pp, process.hrtime.bigint())+"\n======\n";
    console.log(pp_end);
    if (!preprocess_new_result) {
        shell.ShellString(pp_end).toEnd(pre_process_dir + 'logPreProcessOutput');
    }

    return(output_stats);
}

function tremoloIdentification(output_name, output_path, mgf, mz_tol, sim_tol, top_k, verbose, verbose_search) {
    const start_tremolo = process.hrtime.bigint();
    console.log('*** Step 6 - Calling tremolo to perform an in-silico spectral library identification against UNPD *** \n');

    // check if output_path exists
    if (!shell.test('-e', output_path))
    {
        shell.mkdir("-p", output_path);
    }

    // Give the path to the CSV file containing the description of the database
    var db_desc = __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_DB.csv";

    console.log(' Converting the mgf file to a pklbin file \n');
    // run the mgf file converter to pklbin
    const tmpDir = "/tmp/ISDB_tremolo_NP3/" + Math.random().toString(36).substring(2, 10);
    const tremoloTmpDirPath = tmpDir + "/Data/results";
    shell.mkdir("-p", tremoloTmpDirPath);
    var resExec = shell.exec(__dirname+'/src/ISDB_tremolo_NP3/Data/tremolo/convert '+ mgf+
        tremoloTmpDirPath + '/spectra_mgf_'+
        output_name+'.pklbin', {async:false, silent: (verbose <= 0)});

    if (!resExec.code) { // error code is 0
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR. Could not convert the spectra mgf to a pklbin file.\n');
        return(resExec.code);
    } else {
        console.log("DONE!\n");
    }

    shell.ShellString("EXISTING_LIBRARY_MGF="+__dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p01.mgf " +
        __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p02.mgf " +
        __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p03.mgf " +
        __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p04.mgf " +
        __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p05.mgf " +
        __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p06.mgf " +
        __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p07.mgf " +
        __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p08.mgf " +
        __dirname+"/src/ISDB_tremolo_NP3/Data/dbs/UNPD_ISDB_R_p09.mgf\n\n" +
        "searchspectra="+tremoloTmpDirPath+"/spectra_mgf_"+output_name+".pklbin\n\n" +
        "RESULTS_DIR="+tremoloTmpDirPath+"/Results_tremolo_"+output_name+".out\n\n" +
        "tolerance.PM_tolerance="+mz_tol+"\n\n" +
        "search_decoy=0\n\n" +
        "SCORE_THRESHOLD="+sim_tol+"\n\n" +
        "TOP_K_RESULTS="+top_k+"\n\n" +
        "NODEIDX=0\n" +
        "NODECOUNT=1\n\n" +
        "SEARCHABUNDANCE=0\n" +
        "SLGFLOADMODE=1").to(tremoloTmpDirPath+'/scripted_'+output_name+'.params');

    console.log(' Running the tremolo search \n');

    // run the tremolo search
    resExec = shell.exec(__dirname+'/src/ISDB_tremolo_NP3/Data/tremolo/main_execmodule ExecSpectralLibrarySearch ' +
        +tremoloTmpDirPath+'/scripted_'+output_name+'.params ',
        {async:false, silent: (verbose_search === 0)});
    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR. Couldn\'t run the spectral library search.\n');
        return(resExec.code);
    } else {
        console.log("DONE!\n");
        shell.ShellString(resExec.stdout).to(output_path+"/logTremolo")
    }

    resExec = shell.exec(python3()+' '+__dirname+'/src/ISDB_tremolo_NP3/Data/dbs/treat.py ' +
        tremoloTmpDirPath+'/Results_tremolo_' +
        output_name+'.out ' +output_path + ' ' +db_desc,
        {async:false, silent: (verbose <= 0)});
    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log("\nERROR. Couldn't treat the tremolo result file.\n");
    } else {
        console.log("DONE!\n");
    }

    // removing tremolo output temporary directory
    shell.rm("-rf", tmpDir);

    var tremolo_end = "Tremolo search ended!\n======\nFinish Step 6 "+printTimeElapsed_bigint(start_tremolo, process.hrtime.bigint())+"\n======\n";
    console.log(tremolo_end);
    shell.ShellString(tremolo_end).toEnd(output_path+"/logTremolo");
    return(resExec.code);
}

function callCreateBatchLists(metadata, raw_data_path, output_path, output_name, processed_data_name, verbose)
{
    const start_batchlist = process.hrtime.bigint();
    console.log('*** Step 3 - Creating the NP3_MSCluster specification lists for '+output_name+' ***\n');
    try {
        var resExec = shell.exec('Rscript '+__dirname+'/src/create_batch_lists.R ' + metadata + ' ' + raw_data_path + ' ' +
            output_path + ' ' + output_name + ' ' + processed_data_name, {async: false, silent: (verbose <= 0)});
    } catch (e) {
        if (verbose <= 0) {
            console.log(e);
        }
        console.log("\nERROR");
        process.exit(1);
    }
    if (resExec.code) {
        if (verbose <= 0) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
        }
        console.log('\nERROR');
        process.exit(1);
    } else {
        console.log('\nDONE! '+printTimeElapsed_bigint(start_batchlist, process.hrtime.bigint())+"\n");
    }
}

function renameTremoloJoinedIds(path_clean_count, path_tremolo_result, verbose)
{
    const start_renametremolo = process.hrtime.bigint();
    console.log('*** Renaming the tremolo result with the joined IDs ***\n');
    try {
        var resExec = shell.exec('Rscript '+__dirname+'/src/tremolo_update_joinedIds.R ' + path_tremolo_result + ' ' +
            path_clean_count, {async: false, silent: (verbose <= 0)});
    } catch (e) {
        if (verbose <= 0) {
            console.log(e);
        }
        console.log("\nERROR");
    }

    console.log('DONE! '+printTimeElapsed_bigint(start_renametremolo, process.hrtime.bigint())+"\n");
}

function mergeTremoloResults(path_tremolo_result, max_results, path_count_files, verbose)
{
    const start_mergetremolo = process.hrtime.bigint();
    console.log('*** Merging the tremolo result with the counts tables ***\n');
    try {
        shell.exec(python3()+' '+__dirname+'/src/ISDB_tremolo_NP3/Data/dbs/tremolo_merge_count.py ' + path_tremolo_result +
            ' ' + max_results + ' ' + path_count_files.join(' '), {async: false, silent: (verbose <= 0)});
    } catch (e) {
        if (verbose <= 0) {
            console.log(e);
        }
        console.log("\nERROR");
    }

    console.log('\nDONE! '+printTimeElapsed_bigint(start_mergetremolo, process.hrtime.bigint())+"\n");
}

// function callMetfragPubChem(output_name, output_path, method, ion_mode, ppm_tolerance, fragment_tolerance,scale_factor,
//                             verbose)
// {
//     console.log('*** Step 6 - Running a MetFrag identification with PubChem for ' + output_name + ' ***\n');
//
//     resExec = shell.exec('Rscript '+__dirname+'/src/metfrag_autosearch.R ' + output_path + "/" + ' ' +
//         method + ' ' +ion_mode + ' ' + ppm_tolerance + ' ' + fragment_tolerance +
//         ' ' + scale_factor, {async: false, silent: (verbose <= 0)});
//     // { code:..., stdout:... , stderr:... }
//     if (resExec.code) {
//         if (verbose <= 0) {
//             console.log(resExec.stdout);
//             console.log(resExec.stderr);
//         }
//         console.log('ERROR\n');
//     } else {
//         console.log('DONE!\n');
//     }
// }

function checkCountsConsistency(output_path, processed_data, metadata, min_peaks_output,
                                clustering=false, clean=false, merge=false)
{
    output_name = basename(output_path);
    output_path = output_path + '/count_tables/';
    var res_all = "*** checkCountsConsistency of job "+output_name+" ***\n\n";
    // check the clustering counts
    if (clustering)
    {
        console.log('\n*** Testing the consistency of the clustering counts for the job '+output_name+' ***\n');
        resExec = shell.exec('Rscript '+__dirname+'/test/test_counts.R ' + processed_data + ' '+ output_path+output_name+'_spectra.csv ' +
            output_path+output_name+'_peak_area.csv '+metadata +' '+ min_peaks_output, {async: false, silent: false});
        if (resExec.code) {
            console.log('\nERROR');
        }
        res_all = res_all + "*clustering*\n"+resExec.stdout+"\n"+resExec.stderr+"\n"
    }
    if (shell.test("-e", output_path+'clean/') && clean)
    {
        console.log('\n*** Testing the consistency of the clean counts for the job '+output_name+' ***\n');
        resExec  = shell.exec('Rscript '+__dirname+'/test/test_counts.R ' + processed_data + ' '+ output_path+'clean/'+output_name+
            '_spectra_clean_ann.csv '+output_path+'clean/'+output_name+'_peak_area_clean_ann.csv '+
            metadata +' '+ min_peaks_output, {async: false, silent: false});
        if (resExec.code) {
            console.log('\nERROR');
        }
        res_all = res_all + "*clean*\n"+resExec.stdout.toString()+"\n"+resExec.stderr.toString()+"\n"
    }
    if (shell.test("-e", output_path+'merge/') && merge)
    {
        console.log('\n*** Testing the consistency of the merged counts for the job '+output_name+' ***\n');
        resExec = shell.exec('Rscript '+__dirname+'/test/test_counts.R ' + processed_data+' '+ output_path+'merge/'+output_name+
            '_spectra_merged_ann.csv '+output_path+'merge/'+output_name+'_peak_area_merged_ann.csv '+
            metadata +' '+ min_peaks_output, {async: false, silent: false});
        if (resExec.code) {
            console.log('\nERROR');
        }
        res_all = res_all + "*merge*\n"+resExec.stdout+"\n"+resExec.stderr+"\n"
    } else if (merge) {
		console.log('*** Testing the consistency of the merged counts for the job '+output_name+' ***\n');
		console.log('\nDone! No symbolic spectrum was created in merge.\n');
	}

    return res_all;
}

function checkCleanMNConsistency(output_path, sim_tol, mn_tol, rt_tol, mz_tol, top_k, max_component_size,
                                 min_matched_peaks)
{
    output_name = basename(output_path);
    var res_all = "*** checkCleanMNConsistency of job "+output_name+" ***\n\n";

    // check molecular networking consistency
    console.log('\n*** Testing the consistency of the clean and annotated counts and the molecular networking for the job '+output_name+' ***\n');
    resExec = shell.exec('Rscript '+__dirname+'/test/test_clean_mn.R ' + output_path+' '+sim_tol+' '+ mn_tol+' '+ rt_tol+' '+ mz_tol+' '+
        top_k+' '+max_component_size+' '+min_matched_peaks, {async: false, silent: false});
    if (resExec.code) {
        console.log('\nERROR');
    }
    res_all = res_all + resExec.stdout.toString()+
        "\n"+resExec.stderr.toString() + "\n";

    return res_all;
}

function checkMNConsistency(output_path, mn_tol, top_k, max_component_size, min_matched_peaks)
{
    output_name = basename(output_path);
    var res_all = "*** checkMNConsistency of job "+output_name+" ***\n\n";

    // check molecular networking consistency
    console.log('\n*** Testing the consistency of the molecular networking for the job '+output_name+' ***\n');
    resExec = shell.exec('Rscript '+__dirname+'/test/test_mn.R ' + output_path+' '+mn_tol+' '+top_k+' '+
        max_component_size+' '+min_matched_peaks,
        {async: false, silent: false});
    if (resExec.code) {
        console.log('\nERROR');
    }
    res_all = res_all + resExec.stdout.toString()+"\n"+resExec.stderr.toString() + "\n";
    return res_all;
}

function callJoinGNPS(cluster_info_path, result_specnets_DB_path, ms_count_path) {
    // check molecular networking consistency
    console.log('\n*** Joining the GNPS identification to the NP3 counts tables ***\n');
    resExec = shell.exec('Rscript '+__dirname+'/src/join_gnps_identification_result.R \"' + cluster_info_path+'\" '+
        result_specnets_DB_path+' '+ms_count_path, {async: false, silent: false});

    if (resExec.code) {
        console.log('ERROR\n');
    } else {
        console.log('DONE!\n');
    }

    return resExec;
}

function isWindows()
{
    return process.platform === "win32" || process.platform === "win64";
}

function osSep()
{
    if (isWindows())
        return  "\\";
    else
        return "/";
}

function python3()
{
    if (isWindows())
    {
        return 'python'
    } else {
        return 'python3'
    }
}

function defaultModelDir() {
    if (isWindows())
        return __dirname+"\\NP3_MSCluster\\Models_Windows";
    else
        return __dirname+"/NP3_MSCluster/Models";
}

program
    .version('1.1.6',  '--version')
    .usage(' command [options]\n\n' +
        'The NP3 MS workflow is a software system with a collection of scripts to enhance untargeted metabolomics ' +
        'research focused on drug discovery with optimizations towards natural products. \n\n' +
        'The workflow is an automatized procedure to cluster (join) and quantify the MS2 spectra (MS/MS) associated ' +
        'with the same ion, which eluted in concurrent chromatographic peaks (MS1), of a collection of samples ' +
        'from LC-MS/MS experiments. \nIt generates a rank of candidate spectra responsible for the observed hits in ' +
        'bioactivity experiments, suggests the number of metabolites present in the samples and constructs molecular ' +
        'networks to improve the analysis and visualization of the results.\n\n');

program
    .command('setup')
    .description('Check if the NP3 dependencies are installed, try to install missing R and python packages and ' +
        'compile the NP3_MSCluster algorithm\n\n')
    .action(function()
    {
        const start_setup = process.hrtime.bigint();
        console.log('\n*** NP3 workflow setup ***\n\n');
        var resExec;
        var countError = 0;
        var call_cwd = process.cwd();

        console.log('* Merging the UNPD csv file *\n');

        shell.cd(__dirname+'/src/ISDB_tremolo_NP3/Data/dbs');
        resExec = shell.exec('sh merge_db.sh', {async:false});
        if (resExec.code) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
            console.log('\nERROR. Could not merge the UNPD csv file. Tremolo identification is disabled.');
            countError = countError + 1;
        } else {
            console.log("DONE!\n");
        }
        shell.cd(call_cwd);

        // Download the large spec2vec model files - vectors and trainables
        var file_spec2vec = __dirname+'/src/spec2vec_models/spec2vec_UniqueInchikeys_ratio05_filtered_iter_50.model.wv.vectors.npy'
        var url_file_spec2vec = "https://zenodo.org/records/3978054/files/spec2vec_UniqueInchikeys_ratio05_filtered_iter_50.model.wv.vectors.npy?download=1"
        if (!(shell.test('-e', file_spec2vec) && shell.test('-f', file_spec2vec)))
        {
            console.log('* Downloading the spec2vec model vectors table *\n');
            resExec = shell.exec('python '+__dirname+'/src/download_spec2vec_model.py '+
                url_file_spec2vec+' '+file_spec2vec, {async:false, silent: false});
            if (resExec.code) {
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                console.log('\nERROR. Could not download the spec2vec model vectors table.');
                countError = countError + 1;
            } else {
                console.log("DONE!\n");
            }
        }
        file_spec2vec = __dirname+'/src/spec2vec_models/spec2vec_UniqueInchikeys_ratio05_filtered_iter_50.model.trainables.syn1neg.npy'
        url_file_spec2vec = "https://zenodo.org/records/3978054/files/spec2vec_UniqueInchikeys_ratio05_filtered_iter_50.model.trainables.syn1neg.npy?download=1"
        if (!(shell.test('-e', file_spec2vec) && shell.test('-f', file_spec2vec)))
        {
            console.log('* Downloading the spec2vec model trainables table *\n');
            resExec = shell.exec('python '+__dirname+'/src/download_spec2vec_model.py '+
                url_file_spec2vec+' '+file_spec2vec, {async:false, silent: false});
            if (resExec.code) {
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                console.log('\nERROR. Could not download the spec2vec model trainables table.');
                countError = countError + 1;
            } else {
                console.log("DONE!\n");
            }
        }
        
        // check if R is installed
        if (!shell.which('R')) {
            console.log('ERROR. R not found, please ensure R is available and try again.');
            process.exit(1);
        } else {
            // install R packages
            console.log('* Checking R packages requirements *\n');

            resExec = shell.exec('Rscript '+__dirname+'/src/R_requirements.R', {async:false, silent: false});
            if (resExec.code) {
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                console.log('\nERROR. Could not install all R packages, retry with user privileges or do it manually.');

                countError = countError + 1;
            } else {
                console.log("DONE!\n");
            }
        }

        if (isWindows())
        {
            if (!shell.which('python')) {
                console.log('ERROR. Python not found, please ensure python 3 is available and try again.');
                process.exit(1);
            } else {
                if (parseDecimal(shell.exec("python --version", {async:false, silent: true}).stdout.toString().split(' ')[1]) < 3)
                {
                    console.log('ERROR. Python 3 not found, please ensure python 3 is available and try again.');
                    process.exit(1);
                }
            }
        } else {
            if (!shell.which('python3')) {
                console.log('ERROR. Python 3 not found, please ensure python 3 is available and try again.');
                process.exit(1);
            }
        }


        if (!shell.which('pip')) {
            console.log('Pip not found, please ensure pip - the python 3 package management - is available and try again.\n');
            countError = countError + 1;
        } else {
            // install python packages
            console.log('* Checking python 3 packages requirements *\n');
            resExec = shell.exec('pip install -r '+__dirname+'/src/python_requirements --user', {async:false});
            if (resExec.code) {
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                console.log('\nERROR. Could not install all python packages, retry with user privileges or do it manually.');
                countError = countError + 1;
            } else {
                console.log("DONE!\n");
            }
        }

        if (!shell.which('make')) {
            shell.echo('Make not found, please ensure make is available and try again.');
            shell.exit(1);
        }

        // Compile dotproduct with pybind
        console.log('* Compiling the NP3 shifted cosine function for the spectra viewer web app *\n');
        shell.cd(__dirname+'/src/spectra_viewer/src');
        if (!isWindows())
        {
            if (!shell.which('c++')) {
                shell.echo('c++ not found, please ensure C++ is available and try again.');
                shell.exit(1);
            } else {
                resExec = shell.exec('c++ -O3 -Wall -I/usr/include/ -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) norm_dot_product.cpp -o dotprod$(python3-config --extension-suffix)');
                if (resExec.code) {
                    console.log(resExec.stdout);
                    console.log(resExec.stderr);
                    console.log('\nERROR. Could not compile dot product function for spectra viewer web app,\nretry with user privileges or do it manually.\n');
                    countError = countError + 1;
                } else {
                    console.log("DONE!\n");
                }
            }
        } else {
            console.log('Linux/Unix Systems Only, skipped');
        }
        shell.cd(call_cwd);

        console.log('\n* Compiling NP3-MS-Clustering *\n');
        shell.cd(__dirname+'/NP3_MSCluster');

        resExec = shell.exec('make clean', {async:false});
        if (resExec.code) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
            console.log('\nERROR. Could not clean the NP3_MSClustering compilation output. Ensure Make is installed and try again.');
            shell.cd(call_cwd);
            process.exit(1);
        }

        resExec = shell.exec('make', {async:false});
        if (resExec.code) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
            console.log('\nERROR. Could not compile the NP3_MSClustering algorithm. Ensure Make is working and try again.');
            shell.cd(call_cwd);
            process.exit(1);
        } else {
            console.log("DONE!\n\n");
        }

        shell.cd(call_cwd);

        if (countError === 0)
        {
            console.log('NP3 workflow installation complete! ' + printTimeElapsed_bigint(start_setup, process.hrtime.bigint()))
        } else {
            console.log('NP3 workflow installation ended with ' + countError +
                " error(s). Check the error messages, fix the conflicts and retry the setup. " + printTimeElapsed_bigint(start_setup, process.hrtime.bigint()))
        }
    })
    .on('--help', function() {
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow setup');
    });

program
    .command('run')
    .description('Steps 2 to 10: Runs the entire NP3 MS workflow.\n\n')
    .option('-n, --output_name <name>',
        'the job name. It will be used to name the output directory and\n\t\t\t\t\t' +
        'the results from the final clustering integration step. It must have less than 80 characters.\n',
        checkJobNameMaxLength)
    .option('-m, --metadata <file>', 'path to the metadata table CSV file\n')
    .option('-d, --raw_data_path <path>', 'path to the folder containing the input LC-MS/MS raw spectra\n\t\t\t\t\t' +
        'data files (mzXML format is recommended)\n')
    .option('-o, --output_path <path>', 'path to where the output directory will be created\n')
    .option('-f, --fragment_tolerance [x]', 'the tolerance in Daltons for fragment peaks. Peaks in the\n\t\t\t\t\t' +
        'original MS/MS spectra that are closer than this get merged in\n\t\t\t\t\t' +
        'the clustering jobs (Step 3). Also used in the pre process (Step 2), in the\n\t\t\t\t\t' +
        'spectra similarity comparisons and in the cleaning Step 5',
        parseFloat, 0.05)
    .option('-z, --mz_tolerance [x]', 'this is the tolerance in Daltons for the m/z of the\n\t\t\t\t\t' +
        'precursor that determines if two spectra will be compared\n\t\t\t\t\t' +
        'and possibly joined. Used in the clustering jobs (Step 3),\n\t\t\t\t\t' +
        'in the cleaning (Step 5), in the library identifications (Step 6)\n\t\t\t\t\t' +
        'and in the annotation of ionization variants (Step 7)', parseFloat, 0.025)
    .option('-p, --ppm_tolerance [x]', 'the maximal tolerated m/z deviation in parts per million\n\t\t\t\t\t' +
        '(ppm) to be used in the pre processing (Step 2). Typically set to a\n\t\t\t\t\t' +
        'generous multiple of the mass accuracy of the mass\n\t\t\t\t\t' +
        'spectrometer', parseFloat, 15)
    .option('-a, --ion_mode [x]', 'the precursor ion mode. One of the following numeric values\n\t\t\t\t\t' +
        'corresponding to an ion adduct type: \'1\' = [M+H]+ or\n\t\t\t\t\t' +
        '\'2\' = [M-H]-', convertIonMode,1)
    .option('-i, --similarity_function [x]', 'the similarity function to be used in the spectra comparison \n\t\t\t\t\t' +
        'to create the pairwise similarity table after clustering and clean steps. \n\t\t\t\t\t' +
        'If "spec2vec" is selected, the model trained on UniqueInchikey subset (12,797 spectra) \n\t\t\t\t\t' +
        'is used by spec2vec in the spectra comparison and the matchms library is used to compute \n\t\t\t\t\t' +
        'the number of matched peaks between the compared spectra; otherwise, the NP3 shifted cosine \n\t\t\t\t\t' +
        'function is used. One of "np3_shifted_cosine" or "spec2vec".', convertSimilarityFunc, "np3_shifted_cosine")
    .option('-s, --similarity [x]', 'the minimum similarity to be consider in the hierarchical\n\t\t\t\t\t' +
        'clustering of Step 3, starts in 0.70 and decrease to X in 15 rounds.\n\t\t\t\t\t' +
        'Also used in the clean Step 5', parseFloat, 0.55)
    .option('-g, --similarity_blank [x]', 'the minimum similarity to be consider in the hierarchical\n\t\t\t\t\t' +
        'clustering of the blank clustering steps, starts in 0.70\n\t\t\t\t\t' +
        'and decrease to X in 15 rounds. Only used in the \n\t\t\t\t\t' +
        'clustering of blank samples (column SAMPLE\_TYPE equals\n\t\t\t\t\t' +
        '\'blank\' in the metadata table)', parseFloat, 0.3)
    .option('-w, --similarity_mn [x]', 'the minimum similarity score that must occur between a pair\n\t\t\t\t\t' +
        'of consensus spectra to connect them with a link in the\n\t\t\t\t\t' +
        'molecular network of similarity. Lower values will increase the\n\t\t\t\t\t' +
        'components sizes by inducing the connection of less related\n\t\t\t\t\t' +
        'spectra; and higher values will limit the components\n\t\t\t\t\t' +
        'sizes to the opposite', parseFloat, 0.6)
    .option('-t, --rt_tolerance [x,y]', 'tolerances in seconds for the retention time width of the\n\t\t\t\t\t' +
        'precursor that determines if two spectra will be compared\n\t\t\t\t\t' +
        'and possibly joined. It is directly applied to the retention\n\t\t\t\t\t' +
        'time minimum (subtracted) and maximum (added) of the spectra.\n\t\t\t\t\t' +
        'It enlarges the peak boundaries to deal with misaligned samples\n\t\t\t\t\t' +
        'or ionization variant spectra. The first tolerance [x] is used \n\t\t\t\t\t' +
        'in the data, blank and batches integration steps from the clustering\n\t\t\t\t\t' +
        'Step 3 and in Steps 2 (if no previous result is provided) and 7 (chemical annotations); \n\t\t\t\t\t' +
        'and the tolerance [y] is used \n\t\t\t\t\t' +
        'in the final integration step from the clustering Step 3 and in the clean Step 5', splitListFloat, [1,2])
    .option('-k, --net_top_k [x]', 'the maximum number of connections for one single node in the\n\t\t\t\t\t' +
        'molecular network of similarity. A link between two nodes\n\t\t\t\t\t' +
        'is kept only if both nodes are within each other\'s [x]\n\t\t\t\t\t' +
        'most similar nodes. Keeping this value low makes \n\t\t\t\t\t' +
        'very large networks (many nodes) much easier to visualize',15)
    .option('-x, --max_component_size [x]', 'the maximum number of nodes that each component of \n\t\t\t\t\t' +
        'the molecular network of similarity must have (Step 10). The links of \n\t\t\t\t\t' +
        'this network will be removed using an increasing cosine \n\t\t\t\t\t' +
        'threshold until each component has at most X nodes. \n\t\t\t\t\t' +
        'Keeping this value low makes very large networks (many nodes \n\t\t\t\t\t' +
        'and links) much easier to visualize.',200)
    .option('--min_matched_peaks [x]', 'The minimum number of common peaks that two spectra must ' +
        'share to be connected by an edge in the filtered SSMN. Connections ' +
        'between spectra with less common peaks than this cutoff will be ' +
        'removed when filtering the SSMN. Except for when one of the spectra ' +
        'have a number of fragment peaks smaller than the given min_matched_peaks ' +
        'value, in this case the spectra must share at least 2 peaks. ' +
        'The fragment peaks count is performed after the spectra are normalized and cleaned.', parseDecimal, 6)
    .option('--blank_expansion [x]', 'the distance of neighborhood nodes from the blanks in IVAMN to be \n\t\t\t\t\t' +
        'selected for removal in the final protonated networks. \n\t\t\t\t\t' +
        '(0) to only remove blanks nodes,  \n\t\t\t\t\t' +
        '(1) to remove nodes directly connected to a blank node, \n\t\t\t\t\t' +
        '(2 or greater) to remove nodes in a distance equal to 2 or greater from a blank node, \n\t\t\t\t\t' +
        'or (-1) to remove all possible neighbours and ancestors of a blank node (remove blank clusters) from IVAMN ',
        parseDecimal, 0)
    .option('-c, --scale_factor [x]', 'the scaling method to be used in the fragmented peak\'s\n\t\t\t\t\t' +
        'intensities before any dot product comparison (Steps 3 and 5).\n\t\t\t\t\t' +
        'Valid values are: 0 for the natural logarithm (ln)\n\t\t\t\t\t' +
        'of the intensities; 1 for no scaling; and other values\n\t\t\t\t\t' +
        'greater than zero for raising the fragment peaks\n\t\t\t\t\t' +
        'intensities to the power of x (e.g. x = 0.5 is the square\n\t\t\t\t\t' +
        'root scaling). [x] >= 0', parseFloat, 0.5)
    .option('-l, --parallel_cores [x]', 'the number of cores to be used for parallel processing\n\t\t\t\t\t' +
        'in Step 5 spectra comparison. x = 1 for disabling parallelization and x > 2\n\t\t\t\t\t' +
        'for enabling it. x >= 1', parseDecimal, 2)
    .option('-y, --processed_data_name [x]', 'the name of the folder inside the raw_data_path where the\n\t\t\t\t\t' +
        'pre processed data are stored. If the given folder do\n\t\t\t\t\t' +
        'not exists the pre process (Step 2) will be run\n\t\t\t\t\t' +
        'to create it using the default values of the missing Step 2 options.\n\t\t\t\t\t' +
        'Otherwise it will depend on the processed_data_overwrite\n\t\t\t\t\t' +
        'parameter value', "processed_data")
    .option('-q, --processed_data_overwrite [x]', 'A logical "TRUE" or "FALSE" indicating if the pre processed\n\t\t\t\t\t' +
        'data present in the processed_data_name folder (if it\n\t\t\t\t\t' +
        'already exists) should be used (FALSE) or overwritten ' +
        'and pre processed again (TRUE) with the default values of the \n\t\t\t\t\t' +
        'missing Step 2 options', toupper, "FALSE")
    .option('--bflag_cutoff [x]', 'A positive numeric value to scale the interquartile range (IQR)\n\t\t\t\t\t' +
        'of the blank spectra basePeakInt distribution from the clustering result and to allow spectra with a basePeakInt\n\t\t\t\t\t' +
        'value below this distribution median plus IQR*bflag_cutoff to be\n\t\t\t\t\t' +
        'joined with a blank spectrum during the clean Step 5, without relying on the similarity value.\n\t\t\t\t\t' +
        'Or FALSE to disable it.\n\t\t\t\t\t' +
        'The IQR is the range between the 1st quartile (25th quantile) and the\n\t\t\t\t\t' +
        '3rd quartile (75th quantile) of the distribution. The spectra with a\n\t\t\t\t\t' +
        'basePeakInt value <= median + IQR*bflag_cutoff (from the\n\t\t\t\t\t' +
        'blank spectra basePeakInt distribution) and BFLAG TRUE will be joined to a blank spectrum\n\t\t\t\t\t' +
        'in the clean Step 5. This cutoff will affect the spectra with BFLAG TRUE\n\t\t\t\t\t' +
        'that would not get joined to a blank spectra when relying only on the\n\t\t\t\t\t' +
        'similarity cutoff. This is a turn around to the fact that blank spectra\n\t\t\t\t\t' +
        'have low quality spectra and thus can not fully rely on the similarity values.',
        parseBFLAGcutoff, 1.5)
    .option('--noise_cutoff [x]', 'A positive numeric value to scale the interquartile range (IQR)\n\t\t\t\t\t' +
        'of the blank spectra basePeakInt distribution from the clustering Step 3 result and to remove the spectra with a basePeakInt\n\t\t\t\t\t' +
        'value below this distribution median plus IQR*noise_cutoff after the clean Step 5.\n\t\t\t\t\t' +
        'Or FALSE to disable it.\n\t\t\t\t\t' +
        'The IQR is the range between the 1st quartile (25th quantile) and the\n\t\t\t\t\t' +
        '3rd quartile (75th quantile) of the distribution. \n\t\t\t\t\t' +
        'When no blank sample is present in the metadata, the full distribution is used. \n\t\t\t\t\t' +
        'This cutoff will affect the spectra with with a low\n\t\t\t\t\t' +
        'basePeakInt value that probably are noise features. \n\t\t\t\t\t' +
        'If the clustering Step 3 results in more than 25000 spectra, \n\t\t\t\t\t' +
        'the noise cutoff will be applied before the clean Step 5 to prevent a long processing time',
        parseNOISEcutoff, "FALSE")
    .option('-u, --rules [x]', ' path to the CSV file following the NP3 rules table format with the\n\t\t\t\t\t' +
        'accepted ionization modification rules for detecting adducts,\n\t\t\t\t\t' +
        'multiple charge, dimers/trimers and their combination with\n\t\t\t\t\t' +
        'neutral losses. To be used by the annotation algorithm (Step 7)',__dirname+"/rules/np3_modifications.csv")
    .option('-r, --trim_mz [x]', 'A logical "TRUE" or "FALSE" indicating if the spectra fragmented \n\t\t\t\t\t' +
        'peaks around the precursor m/z +-20 Da should be deleted \n\t\t\t\t\t' +
        'before the pairwise comparisons. If "TRUE" this removes the \n\t\t\t\t\t' +
        'residual precursor ion, which is frequently observed in MS/MS \n\t\t\t\t\t' +
        'spectra acquired on qTOFs.', toupper, "TRUE")
    .option('--max_shift [x]', 'Maximum difference between precursor m/zs that will be used in the search of ' +
        'shifted m/z fragment ions in the NP3 shifted cosine function. ' +
        'Shifts greater than this value will be ignored and not used in the cosine computation. ' +
        'It can be useful to deal with local modifications of the same compound.', parseFloat, 200)
    .option('-b, --max_chunk_spectra [x]', "Maximum number of spectra (rows) to be loaded and processed\n\t\t\t\t\t" +
        "in a chunk at the same time. In case of memory issues this\n\t\t\t\t\t" +
        "value should be decreased. To be used in Steps 5, 7 and 10",parseDecimal,3000)
    .option('-e, --method [name]', 'a character string indicating which correlation coefficient\n\t\t\t\t\t' +
        'is to be computed. One of "pearson", "kendall", or\n\t\t\t\t\t' +
        '"spearman" (Step 9)', convertMethodCorr,"spearman")
    // .option('-i, --metfrag_identification [x]', 'a logical "TRUE" or "FALSE" indicating if the MetFrag tool\n\t\t\t\t\t' +
    //     'should be used for spectra identification search against the\n\t\t\t\t\t' +
    //     'PubChem database of the top correlated spectra. \n\t\t\t\t\t' +
    //     'The Metfrag API is not very fast and could raise some errors, \n\t\t\t\t\t' +
    //     'if this happen the user must stop the process (Step 6)', toupper, "FALSE")
    .option('-j, --tremolo_identification [x]', '(not Windows OS\'s) A logical "TRUE" or "FALSE" indicating if\n\t\t\t\t\t' +
        'the Tremolo tool should be used for the spectra matching\n\t\t\t\t\t' +
        'against the ISDB from the UNPD (Step 6)', toupper, "TRUE")
    .option('-v, --verbose [x]', 'for values X>0 show the scripts output information.\n\t\t\t\t\t' +
        ' For values greater or equal to 10 a consistency test of \n\t\t\t\t\t' +
        'the results is also performed.', parseDecimal, 0)
    .action(function(options) {
        //console.log(options);
        // check mandatory params
        if (typeof options.output_name === 'undefined') {
            console.error('\nMissing the mandatory \'output_name\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.metadata === 'undefined') {
            console.error('\nMissing the mandatory \'metadata\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.raw_data_path === 'undefined') {
            console.error('\nMissing the mandatory \'raw_data_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.output_path === 'undefined' || options.output_path === "undefined") {
            console.error('\nMissing the mandatory \'output_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        // check some parameters values
        ["similarity", "similarity_blank", "similarity_mn"].forEach(function(parm)
        {
            var value = options[parm];
            if (isNaN(value) || value > 1 || value < 0)
            {
                console.error('\nERROR. Wrong '+parm+' parameter value. The '+parm+' threshold must be a positive numeric value less than 1.0. Execution aborted.');
                process.exit(1);
            }
        });
        if (isNaN(options.fragment_tolerance) || options.fragment_tolerance < 0)
        {
            console.error('\nERROR. Wrong fragment_tolerance parameter value. The fragment tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        if (isNaN(options.mz_tolerance) || options.mz_tolerance < 0)
        {
            console.error('\nERROR. Wrong mz_tolerance parameter value. The m/z mz_tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        if (isNaN(options.ppm_tolerance) || options.ppm_tolerance < 0)
        {
            console.error('\nERROR. Wrong ppm_tolerance parameter value. The PPM tolerance must be a positive numeric value. Execution aborted.');
            process.exit(1);
        }

        const start_run = process.hrtime.bigint();
        // run workflow
        console.log('*** NP3 MS Workflow RUN for \'' +options.output_name+ '\' - Steps 2 to 10 ***\n');

        // set model dir
        // directory where model files are kept. If running MSCluster on Windows and not from the current
        // directory you should specify the path to 'Models_Windows'
        options.model_dir = defaultModelDir();
        // set num of rounds
        // determines how many rounds are used for the hierarchical clustering. [n] <= 20.
        options.num_rounds = 15;
        // set mixture probability
        // the probability wrongfully adding a spectrum to a cluster. [x] < 0.5. Leave it high because mz and rt must match before similarity comparision
        options.mixture_prob = 0.4;
        // set min_peaks_output
        options.min_peaks_output = 5;

        var output_path = options.output_path+'/'+options.output_name;
        var specs_path = output_path + "/spec_lists";

        // check if the raw data is processed, if not process it with the default parms
        // always call pre process script, let it decide if it needs to perfom some action (missing processed files) or not
        var preprocessing_warning = callPreProcessData(options.output_name, options.metadata,
            options.raw_data_path, options, options.verbose);

        if (shell.test('-e', output_path))
        {
            console.log("RUN ERROR: The provided output path already exists, delete it or change the output directory name and retry. \nPath: " + output_path);
            process.exit();
            //shell.rm('-r', output_path)
        }

        callCreateBatchLists(options.metadata, options.raw_data_path, options.output_path, options.output_name,
            options.processed_data_name, options.verbose);

        // save run parameters values
        shell.ShellString('NP3 MS Workflow - version '+ options.parent.version() +'\n\n').toEnd(output_path+'/logRunParms');
        shell.ShellString('output_name: '+options.output_name + "\n\ncmd: \n\n").toEnd(output_path+'/logRunParms');
        shell.ShellString(options.parent.rawArgs.join(' ')).toEnd(output_path+'/logRunParms');

        // copy rules
        shell.cp(options.rules, output_path);

        callClustering(options, output_path, specs_path);

        // clean output folder
        // remove specs folder
        shell.rm('-rf', specs_path);

        output_path = output_path+"/outs/"+options.output_name;

        // create molecular networking output dir
        shell.mkdir("-p", output_path+ "/molecular_networking/similarity_tables");
        // call pairwise comparison for the clustered spectra
        var out_clustered_spec_comp = callPairwiseComparision(options.output_name, output_path + "/molecular_networking/similarity_tables",
            output_path+"/mgf/", options.fragment_tolerance,
            options.scale_factor, options.trim_mz, options.max_shift,options.parallel_cores,
            options.similarity_function, options.verbose);

        // remove mass dissipation in the clustering from area and spectra count
        resExec = callCleanClusteringCounts(options, output_path, options.mz_tolerance,
            options.rt_tolerance[1], options.fragment_tolerance, '',
            out_clustered_spec_comp);

        var counts_path = output_path+"/count_tables/"+options.output_name;
        //resExec = 0
        if (!resExec) // if the clean was succesful, continue for annotation, correlation and merge
        {
            counts_path = output_path+"/count_tables/clean/"+options.output_name;
            var clean_log_output = output_path+"/count_tables/clean/logCleanOutput";
            // call tremolo with the clean mgf and merge results with clean counts files
            if (!isWindows() && options.tremolo_identification === "TRUE") {
                tremoloIdentification(options.output_name, output_path + "/identifications",
                    output_path+"/mgf/"+options.output_name+"_clean.mgf",
                    options.mz_tolerance,0.2, 10, options.verbose, 0);
                // renameTremoloJoinedIds(counts_path+"_spectra_clean.csv",
                //     output_path+ "/identifications/tremolo_results.csv", options.verbose);
                mergeTremoloResults(output_path + "/identifications/tremolo_results.csv",
                    10, [counts_path+"_spectra_clean.csv",
                        counts_path+"_peak_area_clean.csv"]);
            }
            // annotate spectra variants in the clean counts and create the molecular networking of annotations
            // use the fragment_tolerance in both mz and fragment tolerance and use the default absolute ms2 int cutoff
            resExec = callAnnotateCleanCounts(options, output_path,
                options.mz_tolerance,options.mz_tolerance, options.rt_tolerance[0],
                15);

            if (resExec) {
                // annotation step failed, use the clean tables instead
                // call correlation for the cleaned peak area count and spectra area count
                callComputeCorrelation(options.metadata, counts_path+"_peak_area_clean.csv",
                    options.method, 0, clean_log_output, options.verbose);
                callComputeCorrelation(options.metadata, counts_path+"_spectra_clean.csv",
                    options.method, 0, clean_log_output, options.verbose);
            } else {
                // annotation worked
                // call correlation for the cleaned peak area count and spectra area count
                callComputeCorrelation(options.metadata, counts_path + "_peak_area_clean_ann.csv",
                    options.method, 0, output_path+"/count_tables/clean/logAnnotateOutput",
                    options.verbose);
                callComputeCorrelation(options.metadata, counts_path + "_spectra_clean_ann.csv",
                    options.method, 0, output_path+"/count_tables/clean/logAnnotateOutput",
                    options.verbose);
            }

            callMergeCounts(output_path, options.output_name,
                options.raw_data_path + '/' + options.processed_data_name, options.metadata,
                "TRUE", options.method, options.verbose);
        } else {
            // the clean was not successful, call tremolo for the clustered mgf and corr the not clean counts
            var clustering_log_output = output_path+"/logClusteringOutput";
            // call tremole and merge the results with the clustering counts files
            if (!isWindows() && options.tremolo_identification === "TRUE") {
                // tremolo search for not windows OS
                tremoloIdentification(options.output_name, output_path + "/identifications",
                    output_path+"/mgf/"+options.output_name+"_all.mgf",
                    options.mz_tolerance,0.2, 10, options.verbose, 0);
                mergeTremoloResults(output_path + "/identifications/tremolo_results.csv",
                    10, [counts_path+"_spectra.csv",
                        counts_path+"_peak_area.csv"]);
            }

            // call correlation for the mscluster peak area count
            callComputeCorrelation(options.metadata, counts_path+"_peak_area.csv",
                options.method, 0, clustering_log_output, options.verbose);

            // call correlation for the mscluster spectra count
            callComputeCorrelation(options.metadata, counts_path+"_spectra.csv",
                options.method, 0,  clustering_log_output, options.verbose);
        }
        callCreatMN(output_path, options.similarity_mn, options.net_top_k,
            options.max_component_size, options.min_matched_peaks,
            options.max_chunk_spectra, options.blank_expansion, options.verbose);

        // if (options.metfrag_identification === "TRUE")
        // {
        //     if (options.fragment_tolerance > 0.003) {
        //        console.log('Warning: fragment tolerance is too big for searching on PubMed. Setting it to 0.005. To search with a bigger tolerance value use the separated command.\n');
        //        options.fragment_tolerance = 0.003
        //     }
        //     if (options.ppm_tolerance > 5) {
        //       console.log('Warning: ppm tolerance is too big for searching on PubMed. Setting it to 5. To search with a bigger tolerance value use the separated command.\n');
        //        options.ppm_tolerance = 5
        //     }
        //     callMetfragPubChem(options.output_name, output_path, options.method,
        //         options.ion_mode, options.ppm_tolerance, options.fragment_tolerance,
        //         options.scale_factor, options.verbose);
        // }

        // print the pre-processing warning at the end of the process
        if (preprocessing_warning !== undefined) {
            console.log('*** Pre-processing warning - read the messages below to improve your data processing result ****');
            console.log(preprocessing_warning);
            console.log('*** Done for '+options.output_name+' with one major WARNING in the pre process step. See message above ***');
        } else {
            console.log('*** Done for '+options.output_name+' ***');
        }
        if (resExec) {
            console.log('ERROR in the clean and annotation step, check if this is not wanted.');
        }

        var run_end_msg = "\n\nRUN "+printTimeElapsed_bigint(start_run, process.hrtime.bigint());
        console.log(run_end_msg);
        shell.ShellString(run_end_msg+"\n").toEnd(options.output_path+'/'+options.output_name+'/logRunParms');

        if (options.verbose >= 10) {
            console.log("\n*** TESTING ***\n");

            checkCleanMNConsistency(output_path,options.similarity, options.similarity_mn,
                options.rt_tolerance[1], options.mz_tolerance, options.net_top_k,
                options.max_component_size, options.min_matched_peaks);

            checkCountsConsistency(output_path,
                options.raw_data_path+'/'+options.processed_data_name,
                options.metadata, options.min_peaks_output,
                true, true, true);
        }
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('\n');
        console.log('RESULTS:');
        console.log('\n');
        console.log("A directory inside the *output\_path* named with the *output\_name* containing:\n" +
            "- A copy of the *metadata* and the *rules* files and the command line parameters values used in a file " +
            "named 'logRunParms', for reproducibility\n" +
            "- A folder named 'outs' with the clustering steps results in separate folders containing:\n" +
            "    - A subfolder named 'count_tables' with the Step 4 quantification in CSV tables named as " +
            "'<step\_name>\_(spectra|peak\_area).csv'.\n" +
            "    - Onother subfolder named 'clust' with the clusters membership files (which SCANS or msclusterID were joined)\n" +
            "    - A third subfolder named 'mgf' with the resulting clusters consensus spectra in MGF files\n" +
            "    - A text file named 'logNP3MSClusterOutput' with the NP3\_MSCluster log output.\n\n" +
            "The clustering steps results folders are named as 'B\_\<DATA_COLLECTION_BATCH\>\_\<X\>' where " +
            "*\<DATA\_COLLECTION\_BATCH\>* is the data collection batch number in the metadata file of each group of " +
            "samples and *\<X\>* is 0 if it is the result of a *data clustering step* or 1 if it is the result of a " +
            "*blank clustering step*. The *data collection batch integration step* results are stored in folders named " +
            "as 'B\_\<DATA_COLLECTION_BATCH\>'.\n\n" +
            "The final integration step result is located inside the 'outs' directory in a folder named with the " +
            "*output\_name*. This folder contains the final counts and is where the user will find the final results. " +
            "It also contains the following data:\n" +
            "- Inside the 'count_tables' folder two subfolders named 'clean' and 'merge' containing CSV tables with the " +
            "counts from Steps 5, 7, 8 and 9;\n" +
            "- The 'identifications' folder with the tremolo identification results;\n" +
            "- The 'molecular_networking' folder containing:\n" +
            "    - One subfolder named \"similarity_tables\" with the pairwise similarity tables;\n" +
            "    - Three molecular networks edge files (Steps 7 and 10): \n" +
            "        - the molecular network of annotations named as '\<output\_name\>\_ivamn.selfloops'\n" +
            "        - the complete molecular network of similarity named as '\<output\_name\>\_ssmn\_w\_" +
            "\<similarity_mn\>.selfloops', all links with a similarity value above the cut-off are present\n" +
            "        - the filtered molecular network of similarity named as '\<output\_name\>\_ssmn\_w\_" +
            "\<similarity_mn\>\_k\_\<net\_top\_k\>\_x\_\<max\_component\_size\>.selfloops';\n" +
            "    - One CSV table containing the molecular network of annotations attributes and the assigned [M+H]+ " +
            "representatives, named as " +
            "'<output\_name>_ivamn_attributes.csv'"
        );
        console.log('\n');
        console.log('EXAMPLES:');
        console.log('\n');
        console.log('  $ node np3_workflow run --output_name "test_np3" --output_path "/path/where/the/output/will/be/stored" ' +
            '--metadata "/path/to/the/metadata/file/test_np3_metadata.csv" --raw_data_path "/path/to/the/raw/data/dir" ' +
            '--fragment_tolerance 0.01');
        console.log('');
        console.log('  $ node np3_workflow run -n "test_np3_rt_tol" -o "/path/where/the/output/will/be/stored" ' +
            '-m "/path/to/the/metadata/file/test_np3_metadata.csv" -d "/path/to/the/raw/data/dir" ' +
            '-t 3.5,5 -v 10');
        console.log('');
    });

program
    .command('pre_process')
    .description('Step 2: This command runs the pre-process of the LC-MS/MS raw data. It extracts the list of MS1 peaks ' +
        'with their dimension information (minimum and maximum retention times, the peak area and ID) in each sample, ' +
        'matches the MS2 spectra retention time and precursor m/z against this list and assign to each MS2 spectra a ' +
        'MS1 peak that encompasses it. Additionally, a table with the MS1 peaks without any MS2 spectrum m/z and ' +
        'retention time match are stored in a count table of non-fragmented MS1 peaks.\n\n')
    .option('-n, --data_name <name>', 'the data collection name for verbosity\n')
    .option('-m, --metadata <file>', 'path to the metadata table CSV file\n')
    .option('-d, --raw_data_path <path>', 'path to the folder containing the input LC-MS/MS raw spectra\n\t\t\t\t\t' +
        'data files (mzXML format is recommended)\n')
    .option('-y, --processed_data_name [x]', 'The name of the output directory that should be created\n\t\t\t\t\t' +
        'inside the raw_data_path to store the processed data. When\n\t\t\t\t\t' +
        'not using the default directory, this value must be\n\t\t\t\t\t' +
        'informed in the following *run* or *clustering* commands\n\t\t\t\t\t', "processed_data")
    .option('-t, --rt_tolerance [x]', 'the tolerance in seconds used to enlarge the MS1 peak boundaries\n\t\t\t\t\t' +
        'and accept as a match all MS2 ions that have a retention\n\t\t\t\t\t' +
        'time value that is within a MS1 peak range. This value is\n\t\t\t\t\t' +
        'applied to both sides of the MS1 peaks (RTmin - rt_tolerance\n\t\t\t\t\t' +
        'and RTmax + rt_tolerance). Tries to overcome bad MS1 peak integrations.', parseFloat, 3)
    .option('-z, --mz_tolerance [x]', 'the tolerance in Daltons for matching a MS1 peak m/z \n\t\t\t\t\t'+ 
        'with a MS2 spectrum precursor m/z.', parseFloat, 0.05)
    .option('-p, --ppm_tolerance [x]', 'the maximal tolerated m/z deviation in consecutive MS1 scans in\n\t\t\t\t\t' +
        'parts per million (ppm) for the initial ROI definition of\n\t\t\t\t\t' +
        'the R::xcms::centWave algorithm. Typically set to a\n\t\t\t\t\t' +
        'generous multiple of the mass accuracy of the mass\n\t\t\t\t\t' +
        'spectrometer.', parseFloat, 15)
    .option('-a, --ion_mode [x]', 'the precursor ion mode. One of the following numeric values\n\t\t\t\t\t' +
        'corresponding to a ion adduct type: \'1\' = [M+H]+ or\n\t\t\t\t\t' +
        '\'2\' = [M-H]-', convertIonMode,1)
    .option('-e, --peak_width [X,Y]', 'two numeric values separated by comma without spaces and using\n\t\t\t\t\t' +
        'decimal point equals dot, containing the expected\n\t\t\t\t\t' +
        'approximate peak width in chromatographic space. Given\n\t\t\t\t\t' +
        'as a range (min,max) in seconds. The mean value will be\n\t\t\t\t\t' +
        'used to simulate the width of the fake peaks (see\n\t\t\t\t\t' +
        'documentation)', "2,10")
    .option('-s, --snthresh [x]', 'a numeric defining the signal to noise ratio cutoff. \n\t\t\t\t\t' +
        'A parameter of the R::xcms::centWave algorithm.',
        parseFloat, 0)
    .option('-f, --pre_filter [k,I]', 'two numerics separated by comma "," without spaces and with\n\t\t\t\t\t' +
        'decimal point equals dot "." specifying the prefilter step\n\t\t\t\t\t' +
        'for the first analysis step (MS1 ROI detection) of the R::xcms\n\t\t\t\t\t' +
        'centWave algorithm. Mass traces are only retained if they\n\t\t\t\t\t' +
        'contain at least k peaks with intensity >= I.\n\t\t\t\t\t', "1,750")
    .option('-r, --noise [x]', 'a numeric allowing to set a minimum intensity required for\n\t\t\t\t\t' +
        'MS1 centroids to be considered in the first analysis step\n\t\t\t\t\t' +
        '(centroids with intensity < noise are omitted from ROI\n\t\t\t\t\t' +
        'detection). Used in the R::xcms::centWave algorithm. \n\t\t\t\t\t' +
        'Big values can prevent finding a match between a MS2 spectra \n\t\t\t\t\t' +
        'and a MS1 peak, because less MS1 peaks will be present in the final list', parseFloat, 500)
    .option('-c, --mz_center_fun [name]', 'Name of the function to calculate the m/z center of the\n\t\t\t\t\t' +
        'chromatographic peak in the R::xcms::centWave algorithm. Allowed values are: "wMean":\n\t\t\t\t\t' +
        'intensity weighted mean of the peak\'s m/z values,\n\t\t\t\t\t' +
        '"mean": mean of the peak\'s m/z values, "apex": use\n\t\t\t\t\t' +
        'the m/z value at the peak apex, "wMeanApex3": intensity\n\t\t\t\t\t' +
        'weighted mean of the m/z value at the peak apex and the\n\t\t\t\t\t' +
        'm/z values left and right of it and "meanApex3": mean of\n\t\t\t\t\t' +
        'the m/z value of the peak apex and the m/z values left\n\t\t\t\t\t' +
        'and right of it.', "wMeanApex3")
    .option('-u, --integrate_method [x]', 'Integration method for the R::xcms::centWave algorithm. For\n\t\t\t\t\t' +
        'integrate = 1 peak limits are found through descent on the\n\t\t\t\t\t' +
        'mexican hat filtered data, for integrate = 2 the descent\n\t\t\t\t\t' +
        'is done on the real data. The latter method is more\n\t\t\t\t\t' +
        'accurate but prone to noise, while the former is more\n\t\t\t\t\t' +
        'robust, but less exact.', parseDecimal,2)
    .option('-g, --fit_gauss [x]', "a logical \"TRUE\" or \"FALSE\" indicating whether or not a\n\t\t\t\t\t" +
        "Gaussian should be fitted to each peak in the R::xcms::centWave algorithm.", toupper, "FALSE")
    .option('-j, --max_samples_batch_align [x]', 'The maximum number of not blank samples to be selected\n\t\t\t\t\t' +
        'in the metadata sequence for the batch alignment. Due to memory issues\n\t\t\t\t\t' +
        'this value should not exceed 15 samples to avoid crashing\n\t\t\t\t\t' +
        'the script. If x == 0 disable the alignment.', parseDecimal,3)
    .option('-i, --min_fraction [x]', 'a numeric defining the minimum fraction of samples in at\n\t\t\t\t\t' +
        'least one sample group in which the peaks have to be\n\t\t\t\t\t' +
        'present to be considered as a peak group (feature). Used\n\t\t\t\t\t' +
        'in the alignment process to define hook peaks.\n\t\t\t\t\t', parseFloat, 0.3)
    .option('-w, --bw [x]', 'a numeric defining the bandwidth (standard deviation of the\n\t\t\t\t\t' +
        'smoothing kernel) to be used. This option is passed to\n\t\t\t\t\t' +
        'the density method used in the alignment process.', parseFloat, 2)
    .option('-b, --bin_size [x]', 'a numeric defining the size of the overlapping slices in mz\n\t\t\t\t\t' +
        'dimension. This option is passed to the density method used in the alignment process.', parseFloat, 0.05)
    .option('-x, --max_features [x]', 'a numeric with the maximum number of peak groups to be\n\t\t\t\t\t' +
        'identified in a single mz slice. This option is passed to\n\t\t\t\t\t\n' +
        'the density method used in the alignment process.', parseFloat, 100)
    .option('-q, --processed_data_overwrite [x]', 'A logical "TRUE" or "FALSE" indicating if the pre processed\n\t\t\t\t\t' +
        'data present in the processed_data_name folder should be\n\t\t\t\t\t' +
        'overwritten and pre processed again if it already exists\n\t\t\t\t\t', toupper, "TRUE")
    .option('-v, --verbose [x]', 'for values x>0 show the script output information\n\t\t\t\t\t', parseDecimal,0)
    .action(function(options) {
        // console.log(options)
        if (typeof options.data_name === 'undefined') {
            console.error('\nMissing the mandatory \'data_name\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.metadata === 'undefined') {
            console.error('\nMissing the mandatory \'metadata\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.raw_data_path === 'undefined') {
            console.error('\nMissing the mandatory \'raw_data_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (isNaN(options.rt_tolerance) || options.rt_tolerance < 0)
        {
            console.error('\nERROR. Wrong rt_tolerance parameter value. The retention time tolerance must be a positive numeric value given in seconds. Execution aborted.');
            process.exit(1);
        }
        if (isNaN(options.mz_tolerance) || options.mz_tolerance < 0)
        {
            console.error('\nERROR. Wrong mz_tolerance parameter value. The m/z tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        if (isNaN(options.ppm_tolerance) || options.ppm_tolerance < 0)
        {
            console.error('\nERROR. Wrong ppm_tolerance parameter value. The PPM tolerance must be a positive numeric value. Execution aborted.');
            process.exit(1);
        }

        const start_pp = process.hrtime.bigint();
        // run process peak info and align raw data
        console.log('*** NP3 Pre Process - Step 2 ***\n');

        callPreProcessData(options.data_name, options.metadata, options.raw_data_path,
            options, options.verbose);

        console.log('*** Pre processed samples of '+options.data_name+' ***');
        console.log(printTimeElapsed_bigint(start_pp, process.hrtime.bigint()));
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('DETAILS:');
        console.log('');
        console.log('If *max_samples_batch_align* > 0, align the samples to obtain a good guess for the retention time ' +
            'tolerances between samples of the same data collection batch and between all samples. The alignment do not ' +
            'modify the samples retention time, it is only used to suggest a retention time tolerance value.\n\n' +
            'It generates one MGF by sample containing all the detected MS2 spectra enriched with their respectively ' +
            'matched MS1 peak dimensions, e.g. retention time minimum, maximum, peak area and peak ID (given by the ' +
            'peak detection algorithm). When using the same *raw_data_path* to pre process the LC-MS/MS raw data files ' +
            'using different metadata tables, a different *processed\_data\_name* must be used to avoid overwriting ' +
            'useful files and to store the created MGFs in different folders.');
        console.log('');
        console.log('RESULT:');
        console.log('');
        console.log('A folder named *processed_data_name* inside the *raw_data_path* containing:\n' +
            '- One MGF per sample of the metadata table with the MS2 spectra enriched with their matched MS1 peak ' +
            'chromatographic dimensions. Each pre processed file is named as "\<SAMPLE_CODE\>_peak_info.mgf", where ' +
            '\<SAMPLE_CODE\> is as defined in the metadata file located in the *metadata_path*.\n' +
            '- A count file named "MS1_list_with_MS2.csv" with the quantification by sample of all the MS1 peaks that ' +
            'were assigned to a MS2 ion.\n' +
            '- A count file named "MS1_list_no_MS2.csv" with the quantification by sample of all the MS1 peaks that ' +
            'were not assigned to a MS2 ion, e.g. probably not fragmented MS1 peaks, and thus these peaks will be ' +
            'missing in the final clustering counts. It will help to search for the ion\'s isotopic distributions.\n' +
            '- A log file named "logPreProcessStatisticsWarning" or "logPreProcessStatistics" with the statistics of ' +
            'the rate of MS2 spectra without a MS1 peak correspondence and guidelines to improve this step result.\n' +
            '- A table file named "log_MS2_no_MS1peak_match.csv" with the information of the MS2 spectra that did not' +
            ' have a MS1 peak correspondence.\n' +
            '- If *max_samples_batch_align* > 0, a CSV file named "samples_alignment.csv" is created with the maximum ' +
            'misalignment value in seconds for each data collection batch, as defined in the metadata table, and for ' +
            'all samples together. These values serve as a suggestion for the retention time tolerance to be used in' +
            ' the following steps of the workflow. This alignment process does not change the samples retention time.\n' +
            '- A file named "parameters.csv" with the parameter\'s values used, for reproducibility.'
        );
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('- Pre-processing with the default parameters values and align using at most 6 samples per data collection batch:\n');
        console.log(' $ node np3_workflow.js pre_process \-\-data_name "data_UHPLC_qTOF" \-\-metadata ' +
            '"/path/to/the/metadata/file/test_np3_metadata.csv" \-\-raw_data_path "/path/to/the/raw/data/dir" ' +
            '\-\-max_samples_batch_align 6 \-\-verbose 1');
        console.log('');
        console.log('- Pre-processing without alignment and with a different m/z tolerance:\n');
        console.log(' $ node np3_workflow.js pre_process \-\-data_name "data_UHPLC_qTOF" \-m ' +
            '"/path/to/the/metadata/file/test_np3_metadata.csv" \-d "/path/to/the/raw/data/dir" \-j 0 \-z 0.07');
    });

program
    .command('clustering')
    .description('Steps 3 and 4: This command runs the NP3_MSCluster algorithm to perform the clustering of ' +
        'pre-processed MS/MS data into a collection of consensus spectra. Then, it runs the consensus spectra ' +
        'quantification to count the number of spectra and peak area by sample (clustering counts). If necessary, ' +
        'it runs Step 2. And it can also run the library spectra identifications (Step 6) for the collection of ' +
        'consensus spectra.\n\n')
    .option('-n, --output_name <name>', 'the job name. It will be used to name the output directory and \n\t\t\t\t\t' +
        'the results from the final clustering integration step. It must have less than 80 characters.\n',
        checkJobNameMaxLength)
    .option('-m, --metadata <file>', 'path to the metadata table CSV file\n')
    .option('-d, --raw_data_path <path>', 'path to the folder containing the input LC-MS/MS raw spectra\n\t\t\t\t\t' +
        'data files (mzXML format is recommended)')
    .option('-o, --output_path <path>', 'path to where the output directory will be created\n')
    .option('-f, --fragment_tolerance [x]', 'the tolerance in Daltons for fragment peaks. Peaks in the\n\t\t\t\t\t' +
        'original spectra that are closer than this get merged by\n\t\t\t\t\t' +
        'the NP3_MSCluster algorithm. Also used in the pre process (Step 2)\n', parseFloat, 0.05)
    .option('-z, --mz_tolerance [x]', 'this is the tolerance in Daltons for the m/z of the\n\t\t\t\t\t' +
        'precursor that determines if two spectra will be compared\n\t\t\t\t\t' +
        'and possibly joined. Used in the clustering job and\n\t\t\t\t\t' +
        'in the library identifications (Step 6)', parseFloat, 0.025)
    .option('-p, --ppm_tolerance [x]', 'the maximal tolerated m/z deviation in parts per million (ppm)\n\t\t\t\t\t' +
        'to be used in the pre-processing step if ran\n\t\t\t\t\t', parseFloat, 5)
    .option('-a, --ion_mode [x]', 'the precursor ion mode. One of the following numeric values\n\t\t\t\t\t' +
        'corresponding to a ion adduct type: \'1\' = [M+H]+ or \n\t\t\t\t\t\'2\' = [M-H]-', convertIonMode,1)
    .option('-s, --similarity [x]', 'the minimum similarity to be consider in the hierarchical\n\t\t\t\t\t' +
        'clustering, starts in 0.70 and decrease to X in 15 rounds\n\t\t\t\t\t', parseFloat, 0.55)
    .option('-g, --similarity_blank [x]', 'the minimum similarity to be consider in the hierarchical\n\t\t\t\t\t' +
        'clustering of the blank clustering steps, starts in 0.70\n\t\t\t\t\t' +
        'and decrease to X in 15 rounds. Only used in the \n\t\t\t\t\t' +
        'clustering of blank samples (column SAMPLE\\_TYPE equals\n\t\t\t\t\t' +
        '\'blank\' in the metadata table)', parseFloat, 0.3)
    .option('-t, --rt_tolerance [x,y]', 'tolerances in seconds for the retention time width of the\n\t\t\t\t\t' +
        'precursor that determines if two spectra will be compared\n\t\t\t\t\t' +
        'and possibly joined. It is directly applied to the retention\n\t\t\t\t\t' +
        'time minimum (subtracted) and maximum (added) of the spectra.\n\t\t\t\t\t' +
        'The first tolerance [x] is used in the data, blank and\n\t\t\t\t\t' +
        'batches integration steps of the clustering Step 3 and \n\t\t\t\t\t' +
        'in Step 2 (if no previous result exists or \'processed_data_overwrite\' is TRUE);\n\t\t\t\t\t' +
        'and the tolerance [y] is used in the final integration step\n\t\t\t\t\t' +
        'of the clustering Step 3', splitListFloat, [1,2])
    .option('-y, --processed_data_name [x]', 'the name of the folder inside the raw_data_path where the\n\t\t\t\t\t' +
        'pre processed data are stored. If the given folder do\n\t\t\t\t\t' +
        'not exists the pre process (Step 2) will be run\n\t\t\t\t\t' +
        'to create it using the default values of the missing options.\n\t\t\t\t\t' +
        'Otherwise it will depend on the processed_data_overwrite\n\t\t\t\t\t' +
        'parameter value', "processed_data")
    .option('-q, --processed_data_overwrite [x]', 'A logical "TRUE" or "FALSE" indicating if the pre processed\n\t\t\t\t\t' +
        'data present in the processed_data_name folder (if it\n\t\t\t\t\t' +
        'already exists) should be used (FALSE) or overwritten ' +
        'and pre processed again (TRUE) with the default values of the \n\t\t\t\t\t' +
        'missing Step 2 options', toupper, "FALSE")
    .option('-c, --scale_factor [x]', 'the scaling method to be used in the fragmented peak\'s\n\t\t\t\t\t' +
        'intensities before any dot product comparison (Step 3).\n\t\t\t\t\t' +
        'Valid values are: 0 for the natural logarithm (ln) of the\n\t\t\t\t\t' +
        'intensities; 1 for no scaling; and other values greater\n\t\t\t\t\t' +
        'than zero for raising the fragment peaks intensities to\n\t\t\t\t\t' +
        'the power of x (e.g. x = 0.5 is the square root scaling).\n\t\t\t\t\t' +
        '[x] >= 0', parseFloat, 0.5)
    .option('-e, --method [name]', 'a character string indicating which correlation coefficient is\n\t\t\t\t\t' +
        'to be computed. One of "pearson", "kendall", or "spearman"\n\t\t\t\t\t', convertMethodCorr,"spearman")
    .option('-x, --min_peaks_output [x]', 'the minimum number of fragment peaks that a spectrum must have\n\t\t\t\t\t' +
        'to be outputted after the final clustering step. Spectra\n\t\t\t\t\t' +
        'with less than x fragmented peaks will be discarded. x >= 1\n\t\t\t\t\t',5)
    // .option('-i, --metfrag_identification [x]', 'a logical "TRUE" or "FALSE" indicating if MetFrag tool should\n\t\t\t\t\t' +
    //     'be used for identification search against the PubChem\n\t\t\t\t\t' +
    //     'database of the top correlated spectra.', toupper,"FALSE")
    .option('-j, --tremolo_identification [x]', '(not Windows OS\'s) A logical "TRUE" or "FALSE" indicating if\n\t\t\t\t\t' +
        'the Tremolo tool should be used for the spectral matching\n\t\t\t\t\t' +
        'against the ISDB from the UNPD', toupper, "TRUE")
    .option('-b, --max_chunk_spectra [x]', "Maximum number of spectra to be loaded and processed in a\n\t\t\t\t\t" +
        "chunk at the same time. In case of memory issues this\n\t\t\t\t\t" +
        "value should be decreased",parseDecimal,3000)
    .option('-v, --verbose [x]', 'for values X>0 show the scripts output information\n\t\t\t\t\t', parseDecimal, 0)
    .action(function(options) {
        //console.log(options.output_path)
        // check mandatory params
        if (typeof options.output_name === 'undefined') {
            console.error('\nMissing the mandatory \'output_name\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.metadata === 'undefined') {
            console.error('\nMissing the mandatory \'metadata\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.raw_data_path === 'undefined') {
            console.error('\nMissing the mandatory \'raw_data_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.output_path === 'undefined' || options.output_path === "undefined") {
            console.error('\nMissing the mandatory \'output_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        // check some parameters values
        ["similarity", "similarity_blank"].forEach(function(parm)
        {
            var value = options[parm];
            if (isNaN(value) || value > 1 || value < 0)
            {
                console.error('\nERROR. Wrong '+parm+' parameter value. The '+parm+' threshold must be a positive numeric value less than 1.0. Execution aborted.');
                process.exit(1);
            }
        });
        if (isNaN(options.fragment_tolerance) || options.fragment_tolerance < 0)
        {
            console.error('\nERROR. Wrong fragment_tolerance parameter value. The fragment tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        if (isNaN(options.mz_tolerance) || options.mz_tolerance < 0)
        {
            console.error('\nERROR. Wrong mz_tolerance parameter value. The m/z mz_tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        if (isNaN(options.ppm_tolerance) || options.ppm_tolerance < 0)
        {
            console.error('\nERROR. Wrong ppm_tolerance parameter value. The PPM tolerance must be a positive numeric value. Execution aborted.');
            process.exit(1);
        }

        const start_clust = process.hrtime.bigint();
        // run workflow
        console.log('*** NP3 MS Workflow Clustering for \'' +options.output_name+ '\' - Steps 3 and 4 ***\n');

        // set model dir
        // directory where model files are kept. If running MSCluster on Windows and not from the current
        // directory you should specify the path to 'Models_Windows'
        options.model_dir = defaultModelDir();
        // set num of rounds
        // determines how many rounds are used for the hierarchical clustering. [n] <= 20.
        options.num_rounds = 15;
        // set mixture probability
        // the probability wrongfully adding a spectrum to a cluster. [x] < 0.5
        options.mixture_prob = 0.4;

        var output_path;
        var specs_path;

        if (isWindows())
        {
            output_path = options.output_path+"\\"+options.output_name;
            specs_path = output_path + "\\spec_lists";
        } else {
            output_path = options.output_path+"/"+options.output_name;
            specs_path = output_path + "/spec_lists";
        }

        // check if the raw data is processed, if not process it with the default parms
        // always call pre process script, let it decide if it needs to perfom some action (missing processed files) or not
        var preprocessing_warning = callPreProcessData(options.output_name, options.metadata,
            options.raw_data_path, options, options.verbose);


        if (shell.test('-e', output_path))
        {
            console.log("ERROR: The provided output path already exists, delete it or change the output directory name and retry. \nPath:" + output_path);
            process.exit()
            //shell.rm('-r', output_path)
        }

        callCreateBatchLists(options.metadata, options.raw_data_path, options.output_path, options.output_name,
            options.processed_data_name, options.verbose);

        // save run parameters values
        shell.ShellString('NP3 MS Workflow - version '+ options.parent.version() +'\n\n').toEnd(output_path+'/logRunParms');
        shell.ShellString('output_name: '+options.output_name + "\n\ncmd: \n\n").toEnd(output_path+'/logRunParms');
        shell.ShellString(options.parent.rawArgs.join(' ')).toEnd(output_path+'/logRunParms');

        options.bflag_cutoff = 1.5; // just to plot the vertical lines in the basePeakInt distribution

        // run clustering step
        callClustering(options, output_path, specs_path);

        // tremolo search for not windows OS for the clustered mgf
        if (!isWindows() && options.tremolo_identification === "TRUE")
        {
            tremoloIdentification(options.output_name, output_path + "/outs/" + options.output_name + "/identifications",
                output_path+"/outs/"+options.output_name+"/mgf/"+options.output_name+"_all.mgf",
                options.mz_tolerance,0.2, 10, options.verbose, 0)
            mergeTremoloResults(output_path + "/outs/" + options.output_name + "/identifications/tremolo_results.csv",
                10, [output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_spectra.csv",
                    output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_peak_area.csv"])
        }

        // call correlation for the mscluster peak area count
        callComputeCorrelation(options.metadata, output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_peak_area.csv",
            options.method, 0, output_path+"/outs/"+options.output_name+"/logClusteringOutput",
            options.verbose);

        // call correlation for the mscluster spectra count
        callComputeCorrelation(options.metadata, output_path+"/outs/"+options.output_name+"/count_tables/"+options.output_name+"_spectra.csv",
            options.method, 0,  output_path+"/outs/"+options.output_name+"/logClusteringOutput",
            options.verbose);

        // if (options.metfrag_identification === "TRUE")
        // {
        //     callMetfragPubChem(options.output_name, output_path, options.method,
        //         options.ion_mode, 3, 0.003,
        //         options.scale_factor, options.verbose);
        // }

        // print the pre-processing warning at the end of the process
        if (preprocessing_warning !== undefined) {
            console.log('*** Pre-processing warning - read the messages below to improve your data processing result ****');
            console.log(preprocessing_warning);
            console.log('*** Done for '+options.output_name+' with one major WARNING in the pre-processing step. See message above ***');
        } else {
            console.log('*** Done for '+options.output_name+' ***');
        }

        var clust_end_msg = "\n\nClustering "+printTimeElapsed_bigint(start_clust, process.hrtime.bigint());
        console.log(clust_end_msg);
        shell.ShellString(clust_end_msg+"\n").toEnd(output_path+'/logRunParms');

        if (options.verbose >= 10) {
            console.log("\n*** TESTING ***\n\n");

            checkCountsConsistency(output_path+"/outs/"+options.output_name,
                options.raw_data_path+"/"+options.processed_data_name,
                options.metadata, options.min_peaks_output,
                true, false, false);
        }
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log("A directory inside the *output\_path* named with the *output\_name* containing:\n" +
            "- A copy of the *metadata* file and the command line parameters values used in a file " +
            "named 'logRunParms', for reproducibility\n" +
            "- a folder named 'outs' with the clustering steps results in separate folders containing:\n" +
            "- A subfolder named 'count_tables' with the Step 4 quantifications in CSV tables named as " +
            "'<step\_name>\_(spectra|peak\_area).csv'.\n" +
            "- Onother subfolder named 'clust' with the clusters membership files (which SCANS or msclusterID were joined)\n" +
            "- A third subfolder named 'mgf' with the resulting clusters consensus spectra in MGF files\n" +
            "- A text file named 'logNP3MSClusterOutput' with the NP3\_MSCluster log output.\n\n" +
            "The clustering steps results folders are named as 'B\_\<DATA_COLLECTION_BATCH\>\_\<X\>' where " +
            "*\<DATA\_COLLECTION\_BATCH\>* is the data collection batch number in the metadata file of each group of " +
            "samples and *\<X\>* is 0 if it is a data clustering step or 1 if it is a blank clustering step.\n\n" +
            "The final integration step result is located inside the 'outs' directory in a folder named with the " +
            "*output\_name* and it also contains the tremolo identification results inside the " +
            "'identifications' folder.");
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js clustering --output_name "test_np3" --output_path "/path/where/the/output/will/be/stored" ' +
            '--metadata "/path/to/the/metadata/file/test_np3_metadata.csv" --raw_data_path "/path/to/the/raw/data/dir"');
        console.log('');
        console.log('  $ node np3_workflow.js clustering -n "test_np3_rt_tol" -o "/path/where/the/output/will/be/stored" ' +
            '-m "/path/to/the/metadata/file/test_np3_metadata.csv" -d "/path/to/the/raw/data/dir" ' +
            '-t 3.5,5');
    });


program
    .command('clean')
    .description('Step 5: This command runs the pairwise comparisons of the collection of consensus spectra (if not ' +
        'done yet) and then runs the cleaning of the clustering counts. It also runs Step 7 to annotate possible ion ' +
        'variants using the new clean count tables and to create the molecular network of annotations, and runs Step ' +
        '10 to overwrite any old computation of the molecular network of similarities. It can also run the library ' +
        'spectra identifications (Step 6) for the collection of clean consensus spectra.\n\n')
    .option('-m, --metadata <file>', 'path to the metadata table CSV file')
    .option('-o, --output_path <path>', 'path to the final output data folder, inside the outs directory of the clustering result folder. ' +
        'It should contain the mgf folder, the peak area count CSV and the spectra count CSV. The job name will be extracted from here')
    .option('-y, --processed_data_dir <path>', 'the path to the folder inside the raw data folder where the\n\t\t\t\t\t' +
        'pre processed data (MGFs) are stored.')
    .option('-z, --mz_tolerance [x]', 'the tolerance in Daltons that determines if two spectra will be compared and ' +
        'possibly joined. It is also used in Step 7 to detect possible ion variants.', parseFloat, 0.025)
    .option('-t, --rt_tolerance [x,y]', 'tolerances in seconds for the retention time width of the\n\t\t\t\t\t' +
        'precursor that determines if two spectra will be compared\n\t\t\t\t\t' +
        'and possibly joined. It is directly applied to the retention\n\t\t\t\t\t' +
        'time minimum (subtracted) and maximum (added) of the spectra.\n\t\t\t\t\t' +
        'It enlarges the peak boundaries to deal with disaligned samples\n\t\t\t\t\t' +
        'or ionization variant spectra. The first tolerance [x] is used \n\t\t\t\t\t' +
        'in the annotation Step 7; and the tolerance [y] is used \n\t\t\t\t\t' +
        'in the clean Step 5', splitListFloat, [1,2])
    .option('-a, --ion_mode [x]', 'the precursor ion mode. One of the following numeric values corresponding ' +
        'to a ion adduct type: \'1\' = [M+H]+ or \'2\' = [M-H]-', convertIonMode,1)
    .option('-i, --similarity_function [x]', 'the similarity function to be used in the spectra comparison \n\t\t\t\t\t' +
        'to create the pairwise similarity table after clustering and clean steps. \n\t\t\t\t\t' +
        'If "spec2vec" is selected, the model trained on UniqueInchikey subset (12,797 spectra) \n\t\t\t\t\t' +
        'is used by spec2vec in the spectra comparison and the matchms library is used to compute \n\t\t\t\t\t' +
        'the number of matched peaks between the compared spectra; otherwise, the NP3 shifted cosine \n\t\t\t\t\t' +
        'function is used. One of "np3_shifted_cosine" or "spec2vec".', convertSimilarityFunc, "np3_shifted_cosine")
    .option('-s, --similarity [x]', 'the similarity to be consider when joining clusters and merging ' +
        'their counts in not blank samples', parseFloat, 0.55)
    .option('-g, --similarity_blank [x]', 'the similarity to be consider when joining clusters and ' +
        'merging their counts in blank samples (column SAMPLE\\_TYPE equals \'blank\' in the metadata table)', parseFloat, 0.3)
    .option('-w, --similarity_mn [x]', 'the minimum similarity score that must occur between a pair of ' +
        'consensus MS/MS spectra in order to create an edge in the molecular networking. Lower values will increase the ' +
        'component size of the clusters by inducing the connection of less related MS/MS spectra; and higher values will ' +
        ' limit the components sizes to the opposite', parseFloat, 0.6)
    .option('-f, --fragment_tolerance [x]', 'the tolerance in Daltons for fragment peaks. Peaks in a cluster ' +
        'spectrum that are closer than this are considered the same.', parseFloat, 0.05)
    .option('--bflag_cutoff [x]', 'A positive numeric value to scale the interquartile range (IQR)\n\t\t\t\t\t' +
        'of the blank spectra basePeakInt distribution from the clustering result and to allow spectra with a basePeakInt\n\t\t\t\t\t' +
        'value below this distribution median plus IQR*bflag_cutoff to be\n\t\t\t\t\t' +
        'joined with a blank spectrum during the clean Step 5, without relying on the similarity value.\n\t\t\t\t\t' +
        'Or FALSE to disable it.\n\t\t\t\t\t' +
        'The IQR is the range between the 1st quartile (25th quantile) and the\n\t\t\t\t\t' +
        '3rd quartile (75th quantile) of the distribution. The spectra with a\n\t\t\t\t\t' +
        'basePeakInt value <= median + IQR*bflag_cutoff (from the\n\t\t\t\t\t' +
        'blank spectra basePeakInt distribution) and BFLAG TRUE will be joined to a blank spectrum\n\t\t\t\t\t' +
        'in the clean Step 5. This cutoff will affect the spectra with BFLAG TRUE\n\t\t\t\t\t' +
        'that would not get joined to a blank spectra when relying only on the\n\t\t\t\t\t' +
        'similarity cutoff. This is a turn around to the fact that blank spectra\n\t\t\t\t\t' +
        'have low quality spectra and thus can not fully rely on the similarity values.',
        parseBFLAGcutoff, 1.5)
    .option('--noise_cutoff [x]', 'A positive numeric value to scale the interquartile range (IQR)\n\t\t\t\t\t' +
        'of the blank spectra basePeakInt distribution from the clustering Step 3 result and to remove the spectra with a basePeakInt\n\t\t\t\t\t' +
        'value below this distribution median plus IQR*noise_cutoff after the clean Step 5.\n\t\t\t\t\t' +
        'Or FALSE to disable it.\n\t\t\t\t\t' +
        'The IQR is the range between the 1st quartile (25th quantile) and the\n\t\t\t\t\t' +
        '3rd quartile (75th quantile) of the distribution. \n\t\t\t\t\t' +
        'When no blank sample is present in the metadata, the full distribution is used. \n\t\t\t\t\t' +
        'This cutoff will affect the spectra with with a low\n\t\t\t\t\t' +
        'basePeakInt value that probably are noise features. \n\t\t\t\t\t' +
        'If the clustering Step 3 resulted in more than 25000 spectra, \n\t\t\t\t\t' +
        'the noise cutoff will be applied before the clean Step 5 to prevent a long processing time',
        parseNOISEcutoff, "FALSE")
    .option('-u, --rules [x]', 'path to the CSV file with the accepted ionization modification rules ' +
        'for detecting adducts, multiple charge and dimers/trimers variants, and their combination ' +
        'with neutral losses (Step 7)',__dirname+"/rules/np3_modifications.csv")
    .option('-c, --scale_factor [x]', 'the scaling method to be used in the fragmented peak\'s\n\t\t\t\t\t' +
        'intensities before any dot product comparison (Steps 5 and 7).\n\t\t\t\t\t' +
        'Valid values are: 0 for the natural logarithm (ln) of the\n\t\t\t\t\t' +
        'intensities; 1 for no scaling; and other values greater\n\t\t\t\t\t' +
        'than zero for raising the fragment peaks intensities to\n\t\t\t\t\t' +
        'the power of x (e.g. x = 0.5 is the square root scaling).\n\t\t\t\t\t' +
        '[x] >= 0', parseFloat, 0.5)
    .option('-k, --net_top_k [x]', 'the maximum number of connection for one single node in the\n\t\t\t\t\t' +
        'similarity molecular networking. An edge between two nodes\n\t\t\t\t\t' +
        'is kept only if both nodes are within each other\'s [x]\n\t\t\t\t\t' +
        'most similar nodes. Keeping this value low makes \n\t\t\t\t\t' +
        'very large networks (many nodes) much easier to visualize',15)
    .option('-x, --max_component_size [x]', 'the maximum number of nodes that all component of \n\t\t\t\t\t' +
        'the similarity molecular network must have. The edges of \n\t\t\t\t\t' +
        'this network will be removed using an increasing cosine \n\t\t\t\t\t' +
        'threshold until each network component has at most X nodes. \n\t\t\t\t\t' +
        'Keeping this value low makes very large networks (many nodes \n\t\t\t\t\t' +
        'and edges) much easier to visualize.',200)
    .option('--min_matched_peaks [x]', 'The minimum number of common peaks that two spectra must ' +
        'share to be connected by an edge in the filtered SSMN. Connections ' +
        'between spectra with less common peaks than this cutoff will be ' +
        'removed when filtering the SSMN. Except for when one of the spectra ' +
        'have a number of fragment peaks smaller than the given min_matched_peaks ' +
        'value, in this case the spectra must share at least 2 peaks. ' +
        'The fragment peaks count is performed after the spectra are normalized and cleaned.', parseDecimal, 6)
    .option('--blank_expansion [x]', 'the distance of neighborhood nodes from the blanks in IVAMN to be \n\t\t\t\t\t' +
        'selected for removal in the final protonated networks. \n\t\t\t\t\t' +
        '(0) to only remove blanks nodes,  \n\t\t\t\t\t' +
        '(1) to remove nodes directly connected to a blank node, \n\t\t\t\t\t' +
        '(2 or greater) to remove nodes in a distance equal to 2 or greater from a blank node, \n\t\t\t\t\t' +
        'or (-1) to remove all possible neighbours and ancestors of a blank node (remove blank clusters) from IVAMN ',
        parseDecimal, 0)
    .option('-r, --trim_mz [x]', 'A logical "TRUE" or "FALSE" indicating if the spectra fragmented \n\t\t\t\t\t' +
        'peaks around the precursor m/z +-20 Da should be deleted \n\t\t\t\t\t' +
        'before the pairwise comparisons. If "TRUE" this removes the \n\t\t\t\t\t' +
        'residual precursor ion, which is frequently observed in MS/MS \n\t\t\t\t\t' +
        'spectra acquired on qTOFs.',toupper,"TRUE")
    .option('--max_shift [x]', 'Maximum difference between precursor m/zs that will be used in the search of ' +
        'shifted m/z fragment ions in the NP3 shifted cosine function. ' +
        'Shifts greater than this value will be ignored and not used in the cosine computation. ' +
        'It can be useful to deal with local modifications of the same compound.', parseFloat, 200)
    .option('-l, --parallel_cores [x]', 'the number of cores to be used for parallel processing. ' +
        'x = 1 for disabling parallelization and x > 2 for enabling it. [x] >= 1', parseDecimal, 2)
    .option('-e, --method [name]', 'a character string indicating which correlation coefficient is to be computed. One ' +
        'of "pearson", "kendall", or "spearman"', convertMethodCorr,"spearman")
    // .option('-i, --metfrag_identification [x]', 'a logical "TRUE" or "FALSE" indicating if the MetFrag tool\n\t\t\t\t\t' +
    //     'should be used for identification search against the\n\t\t\t\t\t' +
    //     'PubChem database of the top correlated m/z\'s\n\t\t\t\t\t', toupper,"FALSE")
    .option('-j, --tremolo_identification [x]', 'for not Windows OS\'s. A logical "TRUE" or "FALSE" indicating if the Tremolo tool should be used ' +
        'for the spectral matching against the In-Silico predicted MS/MS spectrum of Natural Products Database (ISDB) from the UNPD ' +
        '(Universal Natural Products Database).', toupper,"TRUE")
    .option('-b, --max_chunk_spectra [x]', "Maximum number of spectra (rows) to be loaded and processed in " +
        "a chunk at the same time. In case of memory issues this value should be decreased",parseDecimal,3000)
    .option('-v, --verbose [x]', 'for values X>0 show the scripts output information', parseDecimal, 0)
    .action(function(options) {
        if (typeof options.metadata === 'undefined') {
            console.error('\nMissing the mandatory \'metadata\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.output_path === 'undefined') {
            console.error('\nMissing the mandatory \'output_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.processed_data_dir === 'undefined') {
            console.error('\nMissing the mandatory \'processed_data_dir\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (isNaN(options.fragment_tolerance) || options.fragment_tolerance < 0)
        {
            console.error('\nERROR. Wrong fragment_tolerance parameter value. The fragment tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        if (isNaN(options.mz_tolerance) || options.mz_tolerance < 0) {
            console.error('\nERROR. Wrong mz_tolerance parameter value. The m/z tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        var output_name = basename(options.output_path);

        // check if pairwise table exists
        var out_clustered_spec_comp = '';
        if (!shell.test('-e', options.output_path + "/molecular_networking/similarity_tables/similarity_table_" +
            output_name + ".csv")) {
            // cread molecular networking output dir
            shell.mkdir("-p", options.output_path + "/molecular_networking/similarity_tables");
            // call pairwise comparison, create folder
            out_clustered_spec_comp = callPairwiseComparision(basename(options.output_path), options.output_path + "/molecular_networking/similarity_tables",
                options.output_path + "/mgf/", options.fragment_tolerance,
                options.scale_factor, options.trim_mz, options.max_shift, options.parallel_cores,
                options.similarity_function, options.verbose);
        }

        const start_clean = process.hrtime.bigint();
        console.log('*** NP3 Clean Clustered Spectra - Step 5 ***\n');

        // remove mass dissipation from the clustering area and spectra count
        // var resExec = callCleanAnnotateClustering(options, options.output_path, options.mz_tolerance,
        //     options.rt_tolerance, options.fragment_tolerance, options.processed_data_dir);
        var resExec = callCleanClusteringCounts(options, options.output_path, options.mz_tolerance,
            options.rt_tolerance[1], options.fragment_tolerance, options.processed_data_dir,
            out_clustered_spec_comp);

        if (!resExec) {   // clean worked, if not windows run tremolo
            var output_clean_path = options.output_path + "/count_tables/clean/";
            // save clean parameters values
            shell.ShellString('output_name: '+output_name + "\n\ncmd: \n\n").toEnd(output_clean_path+'/logCleanParms');
            shell.ShellString(options.parent.rawArgs.join(' ')).toEnd(output_clean_path+'/logCleanParms');

            // run tremolo with the clean mgf and merge results with clean counts files - for not windows OS
            if (!isWindows() && options.tremolo_identification === "TRUE")
            {
                tremoloIdentification(output_name, options.output_path + "/identifications",
                    options.output_path+"/mgf/"+output_name+"_clean.mgf",
                    options.mz_tolerance,0.2, 10, options.verbose, 0);
                // renameTremoloJoinedIds(output_clean_path+ output_name + "_spectra_clean.csv",
                //     options.output_path+ "/identifications/tremolo_results.csv", options.verbose);
                mergeTremoloResults(options.output_path + "/identifications/tremolo_results.csv",
                    10, [output_clean_path + output_name+"_spectra_clean.csv",
                        output_clean_path + output_name+"_peak_area_clean.csv"]);
            }

            // call annotate spectra to identify variant and create the molecular network of annotations
            resExec = callAnnotateCleanCounts(options, options.output_path,
                options.mz_tolerance, options.mz_tolerance, options.rt_tolerance[0],
                15);

            if (resExec) {
                // annotation failed, call correlation in the clean counts
                // call correlation for the cleaned peak area count and spectra area count
                callComputeCorrelation(options.metadata, output_clean_path + output_name + "_peak_area_clean.csv",
                    options.method,0, output_clean_path+"logCleanOutput", options.verbose);

                callComputeCorrelation(options.metadata, output_clean_path + output_name + "_spectra_clean.csv",
                    options.method,0, output_clean_path+"logCleanOutput", options.verbose);
            } else {
                // annotations worked, call correlation in the clean and annotated counts
                // call correlation for the cleaned peak area count and spectra area count
                callComputeCorrelation(options.metadata, output_clean_path + output_name + "_peak_area_clean_ann.csv",
                    options.method,0, output_clean_path+"logAnnotateOutput", options.verbose);

                callComputeCorrelation(options.metadata, output_clean_path + output_name + "_spectra_clean_ann.csv",
                    options.method,0, output_clean_path+"logAnnotateOutput", options.verbose);
            }

            // create MNs
            callCreatMN(options.output_path, options.similarity_mn, options.net_top_k,
                options.max_component_size, options.min_matched_peaks,
                options.max_chunk_spectra, options.blank_expansion, options.verbose);

            // if (options.metfrag_identification === "TRUE")
            // {
            //     callMetfragPubChem(output_name, options.output_path, options.method,
            //         options.ion_mode, 3, options.fragment_tolerance,
            //         options.scale_factor, options.verbose);
            // }

            console.log('*** Done for '+output_name+' ***');
        } else {
            console.log('ERROR in the clean step.');
        }

        console.log("Cleaning "+printTimeElapsed_bigint(start_clean, process.hrtime.bigint()));

        if (options.verbose >= 10 && !resExec ) {
            console.log("\n*** TESTING ***\n\n");

            checkCleanMNConsistency(options.output_path,options.similarity, options.similarity_mn,
                options.rt_tolerance[1], options.mz_tolerance, options.net_top_k,
                options.max_component_size, options.min_matched_peaks);
            // use the default min peaks
            checkCountsConsistency(options.output_path,options.processed_data_dir,
                options.metadata, 5, false, true, false);
        }
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log('One subfolder inside the \'count_tables\' folder is created named \'clean\' containing:\n' +
            '    - A text file named \'analyseCountClusteringClean\' with the clean count analyses\n' +
            '    - Two CSV files with the clustering counts of spectra and peak area cleanned and annotated, named with the ' +
            'suffix \'_clean_ann.csv\'\n' +
            '    - CSV files with the correlation columns added are also included when there is a biocorrelation result\n' +
            'The \'molecular_networking\' folder is also created if not present yet, and inside it: \n' +
            '    - One subfolder named "similarity_tables" with the clean version of the pairwise table of similarity;\n' +
            '    - Two molecular networking\'s edge files: the MN of annotations is named as ' +
            '\'<*output\\_name*>_ivamn.selfloops\' and the MN of similarity is named as ' +
            '\'<*output\\_name*>_ssmn_w_<*similarity_mn*>_k_<*net_top_k*>_x_<*max_component_size*>.selfloop\', where the ' +
            '\'output\\_name\' is extracted from the \'output_path\';\n' +
            '    - One CSV file with the molecular network of annotations edges attributes and protonated representative\n\n'+
            'When running Tremolo the \'identification\' folder is also created (if not present yet) with the' +
            ' identifications results inside it. These identifications are also added as new columns in the c' +
            'reated clean count tables.');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js clean --metadata "/path/to/the/metadata/file/test_np3_metadata.csv" ' +
            '--output_path "/path/to/the/output/dir/test_np3/outs/test_np3" -b 1000');
    });

program
    .command('tremolo')
    .description('Step 6: (for Unix OS and positive ion mode only) This command runs the tremolo tool, used for spectra matching against ' +
        'In-Silico predicted MS/MS spectrum of Natural Products Database (ISDB) from the UNPD (Universal Natural ' +
        'Products Database). It also includes origin and class information of the compounds from NPClassifier, NPAtlas and ClassyFire.\n\n')
    .option('-o, --output_path <path>', 'path to where the spectral library search results will be stored')
    .option('-g, --mgf <path>', 'path to the input MGF file with the MS/MS spectra data to be searched and identified')
    .option('-c, --count_file_path [name]', 'optional paths to the count CSV files, separated by a comma ' +
        'and no space, as outputted by the NP3 workflow where the identifications should be added as new columns. ' +
        'The top_k search results for each msclusterID will be added in 8 identification columns generated by tremolo. ' +
        'The count files header must be in the first row.', splitList, [])
    .option('-z, --mz_tolerance [x]', 'the tolerance for parent mass search in Daltons. Set a small tolerance for ' +
        'desreplication using parent ion mass as prefilter, keeping in mind the resolution of your data. Increase to the ' +
        'wanted range for variable desreplication search (ex: 100 or 200 Da) Caution as this will also increase calculation times!',
        parseFloat, 0.025)
    .option('-s, --similarity [x]', 'the similarity threshold that determines if two spectra are the same. ' +
        'Should be kept low when using in-silico DB, as recommended by the authors. Typically 0.2 to 0.3', parseFloat, 0.2)
    .option('-k, --top_k [x]', 'defines the maximal number of results returned by the tremolo tool', parseDecimal, 10)
    .option('-v, --verbose [x]', 'for values X>0 show the scripts output information', parseDecimal, 0)
    .action(function(options) {
        if (isWindows())
        {
            console.error('Tremolo search is only available for Unix OS.');
            process.exit(1);
        }

        if (typeof options.mgf === 'undefined') {
            console.error('Missing the mandatory \'mgf\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.output_path === 'undefined') {
            console.error('Missing the mandatory \'output_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }

        //call tremoloIdentification(output_name, output_path, count_files, mgf, mz_tol, sim_tol, top_k, verbose, verbose search)
        tremoloIdentification("tremolo_identification", options.output_path, options.mgf,
            options.mz_tolerance, options.similarity, options.top_k, options.verbose,
            options.verbose);
        if (!(options.count_file_path === undefined) && options.count_file_path.length > 0) {
            mergeTremoloResults(options.output_path + "/tremolo_results.csv", options.top_k,
                options.count_file_path);
        }
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log('Two files are created inside the \'output_path\' folder (this folder is created if missing):\n' +
            '    - \'logTremolo\' a text file with the tremolo tool output information;\n' +
            '    - \'tremolo_results.csv\' a CSV file with the tremolo search results for each spectra SCANS present in the provided MGF file.');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js tremolo --output_path "/path/to/the/output/dir/test_np3/outs/test_np3/identifications"' +
            ' --mgf "/path/to/the/output/dir/test_np3/outs/test_np3/mgf/test_np3_all.mgf"');
        console.log('');
        console.log('Concatenate tremolo results in the count files');
        console.log('');
        console.log('  $ node np3_workflow.js tremolo -o "/path/to/the/output/dir/test_np3/outs/test_np3/identifications"' +
            ' -g "/path/to/the/output/dir/test_np3/outs/test_np3/mgf/test_np3_all.mgf" ' +
            '-c "/path/to/the/output/dir/test_np3/outs/test_np3/test_np3_spectra_clean_annotated,' +
            '/path/to/the/output/dir/test_np3/outs/test_np3/test_np3_peak_area_clean_annotated"');
    });

// program
//     .command('metfrag')
//     .description('Step 6: An interactive prompt to run the MetFrag tool for identification search of individual spectra ' +
//         'or of an entire experiment from a MGF file against the PubChem database\n\n')
//     .option('-g, --mgf <path>', 'path to the input MGF file with the MS/MS spectra data to be searched and identified')
//     .option('-o, --output_path <path>', 'path to where the identification results will be saved. ' +
//         'Prefereable inside the final clustering results directory in the \'job_name/outs/job_name/identifications\' folder')
//     .action(function(options) {
//         if (typeof options.mgf === 'undefined') {
//             console.error('Missing the mandatory \'mgf\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
//             process.exit(1);
//         }
//         if (typeof options.output_path === 'undefined') {
//             console.error('Missing the mandatory \'output_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
//             process.exit(1);
//         }
//
//         const { execFileSync } = require('child_process');
//         // run workflow
//         console.log('*** NP3 Comparing Spectra ***');
//
//         try {
//             var resExec = execFileSync("Rscript", ["src/metfrag_interactivesearch.R", options.mgf, options.output_path],
//                 {stdio: 'inherit'});
//         } catch (err) {
//             console.log('\nERROR');
//             //console.log(err.toString().trim());
//             process.exit(1);
//         }
//
//         console.log('\nDONE!\n');
//     })
//     .on('--help', function() {
//         console.log('');
//         console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
//         console.log('');
//         console.log('EXAMPLES:');
//         console.log('');
//         console.log('  $ node np3_workflow.js metfrag --mgf "/path/to/the/mgf/file/input_search.mgf" --output_path "/path/to/the/output/directory/job_name/outs/job_name/identifications"');
//         console.log('');
//         console.log('  $ node np3_workflow.js metfrag --g "/path/to/the/mgf/file/input_search.mgf" -o "/path/to/the/output/directory/job_name/outs/job_name/identifications"');
//     });

program
    .command('annotate_protonated')
    .description('Step 7: (for positive ion mode only) This command runs the annotation of possible ionization variants in the clean count tables ' +
        'and creates the ionization variant annotation molecular network (IVAMN). It searches for adducts, neutral losses, multiple charges, ' +
        'dimers/trimers, isotopes and in-source fragmentation based on numerical equivalences and chemical rules. ' +
        'Finally, it runs a link analysis in the IVAMN to assign some consensus spectra as putative [M+H]+ representatives.\n\n')
    .option('-m, --metadata <file>', 'path to the metadata table CSV file')
    .option('-o, --output_path <path>', 'path to the final output data folder, inside the outs directory of the clustering result folder. ' +
        'It should contain the clean count tables. The job name will be extracted from here')
    .option('-z, --mz_tolerance [x]', 'the tolerance in Daltons for matching the numerical rules of ' +
        'detecting adducts, neutral losses, multiple charge, dimers/trimers and isotopes variants', parseFloat,
        0.025)
    .option('-f, --fragment_tolerance [x]', 'the tolerance in Daltons for matching the numerical rules of ' +
        'in-source fragments and multiple charge isotopic patterns.', parseFloat, 0.025)
    .option('-t, --rt_tolerance [x]', 'tolerance in seconds to enlarge the consensus spectra peak boundaries for detecting ' +
        'concurrent ionization variants spectra.', parseFloat, "1")
    .option('-i, --absolute_ms2_int_cutoff [x]', 'The absolute intensity cutoff for keeping MS2 ' +
        'fragmented peaks with an intensity >= x and to use them in the in-source fragment variant annotation. ' +
        'The MS2 fragmented peaks intensity range from 0 to 1000. (default to 15)', parseFloat, "15")
    .option('-a, --ion_mode [x]', 'the precursor ion mode. One of the following numeric values corresponding ' +
        'to a ion adduct type: \'1\' = [M+H]+ or \'2\' = [M-H]-', convertIonMode,1)
    .option('-u, --rules [x]', 'path to the CSV file with the accepted ionization modification rules for ' +
        'detecting adducts, multiple charge and dimers/trimers variants, and their combination with neutral losses.',
        __dirname+"/rules/np3_modifications.csv")
    .option('-c, --scale_factor [x]', 'the scaling method to be used in the fragmented peak\'s\n\t\t\t\t\t' +
        'intensities before any dot product comparison (Step 5 and 7).\n\t\t\t\t\t' +
        'Valid values are: 0 for the natural logarithm (ln) of the\n\t\t\t\t\t' +
        'intensities; 1 for no scaling; and other values greater\n\t\t\t\t\t' +
        'than zero for raising the fragment peaks intensities to\n\t\t\t\t\t' +
        'the power of x (e.g. x = 0.5 is the square root scaling).\n\t\t\t\t\t' +
        '[x] >= 0', parseFloat, 0.5)
    .option('-b, --max_chunk_spectra [x]', "Maximum number of spectra (rows) to be loaded and processed in " +
        "a chunk at the same time. In case of memory issues this value should be decreased",parseDecimal,3000)
    .option('-v, --verbose [x]', 'for values X>0 show the scripts output information', parseDecimal, 0)
    .action(function(options) {
        if (typeof options.metadata === 'undefined') {
            console.error('\nMissing the mandatory \'metadata\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.output_path === 'undefined') {
            console.error('\nMissing the mandatory \'output_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (isNaN(options.fragment_tolerance) || options.fragment_tolerance < 0)
        {
            console.error('\nERROR. Wrong fragment_tolerance parameter value. The fragment tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        if (isNaN(options.mz_tolerance) || options.mz_tolerance < 0) {
            console.error('\nERROR. Wrong mz_tolerance parameter value. The m/z tolerance must be a positive numeric value given in Daltons. Execution aborted.');
            process.exit(1);
        }
        var output_name = basename(options.output_path);

        // save annotate parameters values
        shell.ShellString('output_name: '+output_name + "\n\ncmd: \n\n").toEnd(options.output_path+'/count_tables/clean/logAnnotateParms');
        shell.ShellString(options.parent.rawArgs.join(' ')).toEnd(options.output_path+'/count_tables/clean/logAnnotateParms');

        const start_ann = process.hrtime.bigint();
        console.log('*** NP3 Annotate Spectra, Create the Molecular Network of Annotations and Assign Protonated Representatives - Step 7 ***\n');

        // call annotate spectra to identify variant and create the molecular network of annotations
        callAnnotateCleanCounts(options, options.output_path,
            options.mz_tolerance, options.fragment_tolerance, options.rt_tolerance,
            options.absolute_ms2_int_cutoff);

        console.log("Annotation [M+H]+ "+ printTimeElapsed_bigint(start_ann, process.hrtime.bigint()));

        if (options.verbose >= 10) {
            console.log("\n*** TESTING ***\n");

            checkMNConsistency(options.output_path, undefined, undefined,
                undefined, undefined);
        }
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log('It creates inside the \'count_tables/clean\' folder:\n ' +
            '- Two CSV files with the clean counts of spectra and peak area concatenated with the annotations ' +
            'result as new columns, named with the suffix \'_annotated.csv\' \n\n' +
            'The \'molecular_networking\' folder is created if not present yet, and inside it is created:\n' +
            '- The molecular network of annotations edge file named as ' +
            '\'\<*output\_name*\>\_ivamn.selfloops\';\n' +
            '- One CSV table containing the molecular network of annotations attributes and the assigned ' +
            '[M+H]+ representatives, named as\'\<*output\_name*\>\_ivamn\_attributes.csv\'\n\n'+
            'Where the \'output\_name\' is extracted from the \'output_path\';');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js annotate_protonated --metadata "/path/to/the/metadata/file/test_np3_metadata.csv" ' +
            '--output_path "/path/to/the/output/dir/test_np3/outs/test_np3" -i 5');
    });

program
    .command('merge')
    .description('Step 8: (for positive ion mode only) This command runs the merge of the clean count tables based on the annotated variants in IVAMN. It ' +
        'creates new symbolic spectra candidates representing the union of each spectra with its annotated variants. ' +
        'By default the merge is only performed for the consensus spectra assigned as a [M+H]+ representative ion, to ' +
        'better account for the quantifications of the true metabolites.\n\n')
    .option('-o, --output_path <path>', 'path to the output data folder, inside the outs directory of the clustering result folder. ' +
        'It should contain the counts_table folder and inside it the clean folder. The job name will be extracted from here')
    .option('-y, --processed_data_dir <path>', 'the path to the folder inside the raw data folder where the\n\t\t\t\t\t' +
        'pre processed data (MGFs) are stored.')
    .option('-m, --metadata <file>', 'path to the metadata table CSV file')
    .option('-p, --merge_protonated [x]', 'A boolean TRUE or FALSE indicating if only the ' +
        '[M+H]+ representative consensus spectra should be merged. If FALSE merge all msclusterID\'s', toupper, "TRUE")
    .option('-e, --method [name]', 'a character string indicating which correlation coefficient is to be computed. One ' +
        'of "pearson", "kendall", or "spearman"', convertMethodCorr,"spearman")
    .option('-v, --verbose [x]', 'for values X>0 show the scripts output information', parseDecimal, 0)
    .action(function(options) {
        if (typeof options.output_path === 'undefined') {
            console.error('\nMissing the mandatory \'output_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.processed_data_dir === 'undefined') {
            console.error('\nMissing the mandatory \'processed_data_dir\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.metadata === 'undefined') {
            console.error('\nMissing the mandatory \'metadata\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }

        // save merge parameters values
        var output_name = basename(options.output_path);
        shell.ShellString('output_name: '+output_name + "\n\ncmd: \n\n").toEnd(options.output_path+'/count_tables/merge/logMergeParms');
        shell.ShellString(options.parent.rawArgs.join(' ')).toEnd(options.output_path+'/count_tables/merge/logMergeParms');

        const start_merge = process.hrtime.bigint();
        // run workflow
        console.log('*** NP3 Merge Counts based on Annotations - Step 8 ***\n');

        // merge double charge and dimer annotation
        callMergeCounts(options.output_path, basename(options.output_path),
            options.processed_data_dir, options.metadata, options.merge_protonated,
            options.method, options.verbose);

        console.log("Merge "+printTimeElapsed_bigint(start_merge, process.hrtime.bigint()));

        if (options.verbose >= 10) {
            console.log("\n*** TESTING ***\n\n");
            // used the default number of min peaks
            checkCountsConsistency(options.output_path,options.processed_data_dir,
                options.metadata, 5, false, false, true);
        }
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log('One subfolder inside the \'count_tables\' folder is created named \'merge\' containing:\n' +
            '- Two CSV files with the cleaned and annotated counts of spectra and peak area merged and new symbolic ' +
            'clusters added as new rows, named with the suffix \'_merged_ann.csv\';\n' +
            '- CSV files with the correlation columns added are also included when there is a biocorrelation result.');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js merge --output_path "/path/to/the/output/dir/test_np3/outs/test_np3" ');
    });

program
    .command('corr')
    .description('Step 9: This command runs the bioactivity correlation to rank the consensus spectra based on the ' +
        'scores computed for the selection of samples and bioactivity values present in the metadata table.\n\n')
    .option('-m, --metadata <file>', 'path to the metadata table CSV file. Used to retrieve the biocorrelation groups')
    .option('-c, --count_file_path <file>', 'path to the count table CSV file')
    .option('-e, --method [name]', 'a character string indicating which correlation coefficient is to be computed. One ' +
        'of "pearson", "kendall", or "spearman"', convertMethodCorr,"spearman")
    .option('-b, --bio_cutoff [name]', 'a bioactivity cutoff value, greater or equal than zero. Bioactivities in the metadata table ' +
        'that are less than the bio_cutoff value will be set to zero.', parseDecimal,0)
    .option('-v, --verbose [x]', 'for values X>0 show the scripts output information', parseDecimal, 0)
    .action(function(options) {
        if (typeof options.metadata === 'undefined') {
            console.error('Missing the mandatory \'metadata\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.count_file_path === 'undefined') {
            console.error('Missing the mandatory \'count_file_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }

        var output_path = basedir(options.count_file_path);
        // save merge parameters values
        shell.ShellString("Corr cmd: \n\n").toEnd(output_path+'/logCorrParms');
        shell.ShellString(options.parent.rawArgs.join(' ')).toEnd(output_path+'/logCorrParms');

        const start_corr = process.hrtime.bigint();
        // run workflow
        console.log('*** NP3 Spectra Count and Bioactivity Correlation - Step 9 ***\n');
        var corr_log_output = output_path+'logCorrOutput';
        // call correlation for the mscluster count
        callComputeCorrelation(options.metadata, options.count_file_path, options.method,
            options.bio_cutoff, corr_log_output, options.verbose);

        console.log("Corr "+ printTimeElapsed_bigint(start_corr, process.hrtime.bigint()));
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log('A CSV table in the count_file directory named as \'\<count_file\>_corr_\<method\>.csv\', ' +
            'without the count_file extension. The new tables contain the original count tables plus for each correlation ' +
            'group and bioactivity score, as defined in the metadata table and new columns with the correlation scores of ' +
            'each consensus spectrum. Also, in a copy of this table, named with the suffix \'_bioAct.csv\', new rows are added with the bioactivities scores of each sample placed above the original ' +
            'counts table header (first rows).\n' +
            '\n' +
            'In the metadata table the columns with a bioactivity score must be named with the prefix "BIOACTIVITY_" and ' +
            'the columns with the correlation groups (samples to be used in each correlation) must be named with the prefix "COR_".\n' +
            'The correlation score will produce NA values with warnings if the counts of the selected samples are all equal 0 or ' +
            'if the bioactivity of the selected samples are all the same, and it will produce "CTE" if the counts of the selected samples ' +
            'have the same values (standard deviation equals 0). ');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js corr --metadata "/path/to/the/metadata/file/test_np3_metadata.csv" --count_file "/path/to/the/metadata/file/np3_job_spectra.csv"')
    });

program
    .command('mn')
    .description('Step 10: This command runs the creation of a spectra similarity molecular network (SSMN) ' +
        'which connects spectra based on the pairwise spectra similarity value above the given similarity cut-off. ' +
        'Then, a filter is applied on this network to remove links between spectra that have less peaks in common than the minimum number of matched peaks, ' +
        'to limit the number of neighbors of each node (number of links) to the top K most similar ones and ' +
        'to limit the size of the components to a maximum number of nodes. The final filtered SSMN contains components that ' +
        'represent the most analogous spectra, possible connecting spectra from similar chemical classes. ' +
        'At the end, a [M+H]⁺ analysis is executed if IVAMN is present, and results in the protonated IVAMN and SSMN. \n\n')
    .option('-o, --output_path <path>', 'path to the output data folder, inside the outs directory of the clustering result folder. ' +
        'It should contain the molecular_networking folder and inside it the similarity_tables folder. The job name will be extracted from here')
    .option('-w, --similarity_mn [x]', 'the minimum similarity score that must occur between a pair of ' +
        'consensus MS/MS spectra in order to create an edge in the molecular networking. Lower values will increase the ' +
        'component size of the clusters by inducing the connection of less related MS/MS spectra; and higher values will ' +
        ' limit the components sizes to the opposite', parseFloat, 0.6)
    .option('-k, --net_top_k [x]', 'the maximum number of connection for one single node in the\n\t\t\t\t\t' +
        'similarity molecular networking. An edge between two nodes\n\t\t\t\t\t' +
        'is kept only if both nodes are within each other\'s [x]\n\t\t\t\t\t' +
        'most similar nodes. Keeping this value low makes \n\t\t\t\t\t' +
        'very large networks (many nodes) much easier to visualize',15)
    .option('-x, --max_component_size [x]', 'the maximum number of nodes that all component of \n\t\t\t\t\t' +
        'the similarity molecular network must have. The edges of \n\t\t\t\t\t' +
        'this network will be removed using an increasing cosine \n\t\t\t\t\t' +
        'threshold until each network component has at most X nodes. \n\t\t\t\t\t' +
        'Keeping this value low makes very large networks (many nodes \n\t\t\t\t\t' +
        'and edges) much easier to visualize.', 200)
    .option('--min_matched_peaks [x]', 'The minimum number of common peaks that two spectra must ' +
        'share to be connected by an edge in the filtered SSMN. Connections ' +
        'between spectra with less common peaks than this cutoff will be ' +
        'removed when filtering the SSMN. Except for when one of the spectra ' +
        'have a number of fragment peaks smaller than the given min_matched_peaks ' +
        'value, in this case the spectra must share at least 2 peaks. ' +
        'The fragment peaks count is performed after the spectra are normalized and cleaned.', parseDecimal, 6)
    .option('--blank_expansion [x]', 'the distance of neighborhood nodes from the blanks in IVAMN to be \n\t\t\t\t\t' +
        'selected for removal in the final protonated networks. \n\t\t\t\t\t' +
        '(0) to only remove blanks nodes,  \n\t\t\t\t\t' +
        '(1) to remove nodes directly connected to a blank node, \n\t\t\t\t\t' +
        '(2 or greater) to remove nodes in a distance equal to 2 or greater from a blank node, \n\t\t\t\t\t' +
        'or (-1) to remove all possible neighbours and ancestors of a blank node (remove blank clusters) from IVAMN ',
        parseDecimal, 0)
    .option('-b, --max_chunk_spectra [x]', "Maximum number of spectra (rows) to be loaded and processed in " +
        "a chunk at the same time. In case of memory issues this value should be decreased",parseDecimal,3000)
    .option('-v, --verbose [x]', 'for values X>0 show the scripts output information', parseDecimal, 0)
    .action(function(options) {
        if (typeof options.output_path === 'undefined') {
            console.error('\nMissing the mandatory \'output_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }

        var output_name = basename(options.output_path);
        // save mn parameters values
        shell.ShellString('output_name: '+output_name + "\n\ncmd: \n\n").toEnd(options.output_path+'/molecular_networking/logMnParms');
        shell.ShellString(options.parent.rawArgs.join(' ')).toEnd(options.output_path+'/molecular_networking/logMnParms');

        const start_mn = process.hrtime.bigint();
        // run workflow
        console.log('*** NP3 Molecular Networking Creation - Step 10 ***\n');
        callCreatMN(options.output_path, options.similarity_mn,
            options.net_top_k, options.max_component_size, options.min_matched_peaks,
            options.max_chunk_spectra, options.blank_expansion, options.verbose);

        console.log("MN "+printTimeElapsed_bigint(start_mn, process.hrtime.bigint()));

        if (options.verbose >= 10) {
            console.log("\n*** TESTING ***\n\n");

            checkMNConsistency(options.output_path, options.similarity_mn, options.net_top_k,
                options.max_component_size, options.min_matched_peaks);
        }
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log('Two files with the molecular networks of similarity are created inside the \'molecular_networking\' folder:\n' +
            '- One named \'\<*output\_name*\>\_ssmn\_w\_\<*similarity_mn*\>.selfloops\' with the complete network, all links with a similarity value above the cut-off are present\n' +
            '- Another named \'\<*output\_name*\>\_ssmn\_w\_\<*similarity_mn*\>\_k\_' +
            '\<*net\_top\_k*\>\_x\_\<*max\_component\_size*\>.selfloops\' with the filtered network;\n' +
            '    \n' +
            'Where the \'output\_name\' is extracted from the \'output_path\';');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js mn --output_path "/path/to/the/output/dir/test_np3/outs/test_np3" ' +
            '--similarity_mn 0.7 --verbose 1');
    });

program
    .command('gnps_result')
    .description('This command join the GNPS library identification result from the Molecular Networking (download clustered spectra) or ' +
        'the Library Search (download all identifications) workflows to the count tables of the NP3 clustering or clean steps\n\n')
    .option('-i, --cluster_info_path [path]', 'If joining the result of a Molecular Networking job, ' +
        'this should be the path to the file inside the folder named ' +
        '\'clusterinfo\' of the downloaded output from GNPS. Not used for results coming from the Library Search workflow.', "")
    .option('-s, --result_specnets_DB_path <path>', 'If joining the result of a Molecular Networking job,' +
        ' this should be the path to the file inside the folder named ' +
        '\'result_specnets_DB\'; and if this is the result of a Library Search workflow,' +
        'this should be the path to the file inside the downloaded folder')
    .option('-c, --count_file_path <path>', 'Path to any of the count tables (peak_area or spectra) resulting ' +
        'from the NP3 clustering or clean steps. If the peak_area is informed and the spectra table file '+
        'exists in the same path (or the opposite), it will merge the GNPS results to both files')
    .action(function(options) {
        if (typeof options.cluster_info_path === 'undefined') {
            console.error('\nMissing the mandatory \'cluster_info_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.result_specnets_DB_path === 'undefined') {
            console.error('\nMissing the mandatory \'result_specnets_DB_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.count_file_path === 'undefined') {
            console.error('\nMissing the mandatory \'count_file_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }

        const start_gnpsjoin = process.hrtime.bigint();
        // run workflow
        console.log('*** Join of the GNPS identification result to the NP3 count files ***\n');

        callJoinGNPS(options.cluster_info_path, options.result_specnets_DB_path,
            options.count_file_path);

        console.log("GNPS_result "+printTimeElapsed_bigint(start_gnpsjoin, process.hrtime.bigint()));
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log('The following columns with the GNPS results are added to the count tables: "gnps_SpectrumID", ' +
            '"gnps_Adduct", "gnps_Smiles", "gnps_CAS_Number", "gnps_Compound_Name", "gnps_LibMZ", "gnps_MZErrorPPM", ' +
            '"gnps_MQScore", "gnps_LibraryQualityString", "gnps_SharedPeaks", "gnps_Organism", "gnps_superclass", ' +
            '"gnps_class", "gnps_subclass" and "gnps_Ion_Source_Instrument". See the GNPS documentation for the ' +
            'description of these columns.\n' +
            '\n' +
            'If there is more than one GNPS result for a single msclusterID the results are concatenated with a ' +
            '\';\', except for the "gnps_Smiles" column which is concatenated with a \',\'.');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js gnps_result --cluster_info_path "/path/to/the/output/dir/GNPS_result/clusterinfo/file" ' +
            '--result_specnets_DB_path "/path/to/the/output/dir/GNPS_result/result_specnets_DB/file.tsv" ' +
            '--ms_count_path "/path/to/the/output/NP3/count_files/count.csv"');
    });

program
    .command('chr')
    .description('This command runs an interactive prompt to extract chromatogram(s) from raw MS1 data files (mzXML, ' +
        'mzData and mzML) and to save to PNG image files. Depending on the provided parameters this can be a total ' +
        'ion chromatogram (TIC - default), a base peak chromatogram (BPC) or an extracted ion chromatogram (XIC) ' +
        'extracted from each sample/file.\n\n')
    .option('-n, --data_name [name]', 'the data collection name. Used for verbosity', "X")
    .option('-m, --metadata <file>', 'path to the metadata table CSV file')
    .option('-d, --raw_data_path <path>', 'path to the folder containing the input LC-MS/MS raw spectra ' +
        'data in mzXML, mzData and mzML format')
    .action(function(options) {
        if (typeof options.metadata === 'undefined') {
            console.error('Missing the mandatory \'metadata\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }
        if (typeof options.raw_data_path === 'undefined') {
            console.error('Missing the mandatory \'raw_data_path\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
            process.exit(1);
        }

        const { execFileSync } = require('child_process');
        // run workflow
        console.log('*** NP3 Extracting Chromatograms ***');

        // call extract chromatogram using xcms
        try {
            var resExec = execFileSync("Rscript",
                [__dirname+"/src/chromatogram_xcms.R", options.data_name, options.metadata,
                options.raw_data_path], {stdio: 'inherit'});
        } catch (err) {
            console.log('\nERROR');
            //console.log(err.toString().trim());
            process.exit(1);
        }

        console.log('\nDONE!\n');
    })
    .on('--help', function() {
        console.log('');
        console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
        console.log('');
        console.log('RESULTS:');
        console.log('');
        console.log('PNG image(s) file(s) with the extracted chromatogram(s) of the selected sample(s) and grouped as specified in the options.');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ node np3_workflow.js chr --metadata "/path/to/the/metadata/file/test_np3_metadata.csv" --data_name "data_np3" --raw_data_path "/path/to/the/raw/data/directory"');
        console.log('');
        console.log('  $ node np3_workflow.js chr -m "/path/to/the/metadata/file/test_np3_metadata.csv" -d "/path/to/the/raw/data/directory"');
    });

// program
//     .command('compare_spectra')
//     .description('An interactive prompt to compare two spectra from a MGF file, to plot them against it other and to ' +
//         'save the image to a PNG file.\n\n')
//     .option('-n, --data_name [name]', 'the data collection name. Used for verbosity', "-")
//     .option('-g, --mgf <path>', 'path to the input MGF file with the MS/MS spectra data to be compared')
//     .action(function(options) {
//         if (typeof options.mgf === 'undefined') {
//             console.error('Missing the mandatory \'mgf\' parameter. See --help for the list of mandatory parameters indicated by angled brackets (e.g. <value>).');
//             process.exit(1);
//         }
//
//         const { execFileSync } = require('child_process');
//         // run workflow
//         console.log('*** NP3 Comparing Spectra ***');
//
//         try {
//             var resExec = execFileSync("Rscript", ["src/compare_spectra.R", options.data_name, options.mgf], {stdio: 'inherit'});
//         } catch (err) {
//             console.log('\nERROR');
//             //console.log(err.toString().trim());
//             process.exit(1);
//         }
//
//         console.log('\nDONE!\n');
//     })
//     .on('--help', function() {
//         console.log('');
//         console.log('Angled brackets (e.g. <x>) indicate required input. Square brackets (e.g. [y]) indicate optional input.');
//         console.log('');
//         console.log('EXAMPLES:');
//         console.log('');
//         console.log('  $ node np3_workflow.js chr --metadata "/path/to/the/metadata/file/test_np3_metadata.csv" --data_name "data_np3" --raw_data_path "/path/to/the/raw/data/directory"');
//         console.log('');
//         console.log('  $ node np3_workflow.js chr -m "/path/to/the/metadata/file/test_np3_metadata.csv" -d "/path/to/the/raw/data/directory"');
//     });

program
    .command('spectra_viewer')
    .description('(for Unix OS only) This command runs an interactive Web App to visualize and compare MS2 spectra. ' +
        'It receives as input a MGF file or a peak list. It is also possible to manipulate, filter, calculate ' +
        'similarity of the spectra and save PNG or SVG plots.\n\n')
    .option('-p, --port [port_number]', 'localhost server port number', "8501")
    .action(function(options) {
        var call_cwd = process.cwd();
        // run workflow
        console.log('*** NP3 Spectra Viewer (press control + c to cancel the operation) ***');
        shell.cd(__dirname+'/src/spectra_viewer');
        try {
            var resExec = shell.exec('streamlit run main.py --server.port '+options.port,{async:false, silent:false});
        } catch (err) {
            console.log(resExec.stdout);
            console.log(resExec.stderr);
            console.log('\nERROR');
            process.exit(1);
        }
        shell.cd(call_cwd);
        console.log('\nDONE!\n');
    })
    .on('--help', function() {
        console.log('');
        console.log('Spectra viewer is a web app running at streamlit framework');
        console.log('');
        console.log('Use -p parameter to change the localhost port');
        console.log('');
        console.log('EXAMPLES:');
        console.log('');
        console.log('  $ spectra_viewer');
        console.log('');
        console.log('  $ spectra_viewer -p 8080');
    });

// TODO: add trycatch in each test case
program
    .command('test')
    .description('This command runs some use cases to test the $NP^{3}$ MS workflow consistency in all steps. ' +
        'This option is intended for debugging purposes, and is not a part of the analysis workflow.\n\n')
    .option('-p, --pre_process [x]', '\'TRUE\' or \'FALSE\' to test the pre process step.', toupper, "FALSE")
    .option('-t, --tremolo [x]', '\'TRUE\' or \'FALSE\' to test the tremolo step. If on Windows OS this should be FALSE.', toupper, "TRUE")
    .option('-s, --skip [x]', 'Skip to test x', 1)
    .action(function(options) {
        const start_test = process.hrtime.bigint();
        var np3_js_call = 'node '+ __dirname +'/np3_workflow.js';

        var unit_test_res = ['@@@@@ UNIT TEST NP3 Shifted cosine @@@@@\n'];
        var test_res = ['@@@@@ TEST 1 @@@@@\n','@@@@@ TEST 2 @@@@@\n','@@@@@ TEST 3 @@@@@\n','@@@@@ TEST 4 @@@@@\n','@@@@@ TEST 5 @@@@@',
            '@@@@@ TEST 6 @@@@@','@@@@@ TEST 7 @@@@@','@@@@@ TEST 8 @@@@@','@@@@@ TEST 9 @@@@@',
            '@@@@@ TEST 10 @@@@@'];
        var resExec;

        // start with the unit tests
        console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
        console.log("@@@@@@ Unit Test - NP3 Shifted Cosine @@@@@");
        console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        resExec = shell.exec('Rscript '+__dirname+'/test/test_np3_shifted_cosine.R',
            {async:false, silent:true});
        if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
            // in case of error show all the emitted msgs
            console.log(resExec.stdout);
            console.log(resExec.stderr);
            unit_test_res[0] = unit_test_res[0] + '\n\nEXEC ERROR';
            console.log('ERROR\n');
        } else {
            unit_test_res[0] = unit_test_res[0] + resExec.stdout;
            console.log('DONE!\n');
        }


        // # not using parallel in the pairwise
        if (options.skip <= 1) {
            shell.rm('-rf', __dirname+'/test/L754_bacs/L754_bacs_all');
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@@ Test 1 - L754_bacs_all - run @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' run -n L754_bacs_all ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata.csv ' +
                '-d '+__dirname+'/test/L754_bacs/mzxml -o '+__dirname+'/test/L754_bacs ' +
                '-j '+options.tremolo+' -v 10 -q '+options.pre_process+' -l 1 '+
                '--noise_cutoff FALSE',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[0] = test_res[0] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[0] = test_res[0] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }
            //console.log("\n\n");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@@ Test 1.1 - L754_bacs_all - gnps_result - Molecular Networking @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' gnps_result ' +
                '-i '+__dirname+'/test/L754_bacs/ProteoSAFe-METABOLOMICS-SNETS-MOLECULARNETWORKING-V2-2dfe22ff-download_clustered_spectra/clusterinfo/0e83d32ce4414494ad9cc12ad3db4824.clusterinfo ' +
                '-s '+__dirname+'/test/L754_bacs/ProteoSAFe-METABOLOMICS-SNETS-MOLECULARNETWORKING-V2-2dfe22ff-download_clustered_spectra/result_specnets_DB/31ba0709274e450295c6da492030f356.tsv ' +
                '-c '+__dirname+'/test/L754_bacs/L754_bacs_all/outs/L754_bacs_all/count_tables/clean/L754_bacs_all_peak_area_clean_ann.csv',
                {async:false, silent:true});
            var gnps_result_mn = "";
            if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emitted msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[0] = test_res[0] + '\n\ngnps_result - Molecular Networking - EXEC ERROR';
                gnps_result_mn = "ERROR";
                console.log('ERROR\n');
            } else {
                gnps_result_mn = resExec.stdout.split('Time Elapsed:')[0];
                // test_res[0] = test_res[0] + '\n*** TESTING - gnps_result - Molecular Networking ***\n\n' + resExec.stdout;
                console.log('DONE!\n');
            }
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@@ Test 1.2 - L754_bacs_all - gnps_result - Library Search @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' gnps_result ' +
                '-s '+__dirname+'/test/L754_bacs/ProteoSAFe-MOLECULAR-LIBRARYSEARCH-V2-da67f38d-download_all_identifications/MOLECULAR-LIBRARYSEARCH-V2-da67f38d-download_all_identifications-main.tsv ' +
                '-c '+__dirname+'/test/L754_bacs/L754_bacs_all/outs/L754_bacs_all/count_tables/clean/L754_bacs_all_spectra_clean_ann.csv',
                {async:false, silent:true});
            var gnps_result_ls = "";
            if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[0] = test_res[0] + '\n\ngnps_result - Library Search - EXEC ERROR';
                gnps_result_ls = "ERROR";
                console.log('ERROR\n');
            } else {
                gnps_result_ls = resExec.stdout.split('Time Elapsed:')[0];
                // test_res[0] = test_res[0] + '\n*** TESTING - gnps_result - Library Search ***\n\n' + resExec.stdout;
                console.log('DONE!\n');
            }
            if (gnps_result_mn == gnps_result_ls && gnps_result_mn != "ERROR" && gnps_result_ls != "ERROR") {
                test_res[0] = test_res[0] +
                    '\n*** Testing - gnps_result - Library Search & Molecular Networking ***\n\nEquality OK!\nDONE! :)\n'
            }
        }

        // # test metadata with more than 10 samples in more than 10 data collections with all types
        // use the same pre processed data
        if (options.skip <= 2) {
            shell.rm('-rf', __dirname+'/test/L754_bacs/L754_bacs_multi_collection_11');
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@@ Test 2 - L754_bacs_multi_collection_11 - run @@@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec=shell.exec(np3_js_call+' run -n L754_bacs_multi_collection_11 ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata_multi_collection_11.csv ' +
                '-d '+__dirname+'/test/L754_bacs/mzxml -o '+__dirname+'/test/L754_bacs/ -j '+options.tremolo+' -v 11 ' +
                '--noise_cutoff 2',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[1] = test_res[1] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[1] = test_res[1] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }
            //console.log("\n\n");
        }

        // # test metadata with all samples in one data collection batch and without blanks
        // # split the run cmd in the sub cmds calls and using a smaller chunk size
        if (options.skip <= 3) {
            shell.rm('-rf', __dirname+'/test/L754_bacs/L754_bacs_one_collection');
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 3 - L754_bacs_one_collection - spec2vec - run @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' run -n L754_bacs_one_collection ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata_one_collection.csv ' +
                '-d '+__dirname+'/test/L754_bacs/mzxml ' +
                '-o '+__dirname+'/test/L754_bacs/ -y processed_data_one_collection -j '+options.tremolo+' -b 200 ' +
                '-v 11 -t 5,10 -q '+options.pre_process+
                ' --bflag_cutoff 1.5 --noise_cutoff 1.5 -i spec2vec',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[2] = test_res[2] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[2] = test_res[2] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }
            //console.log("\n\n");
        }

        // # test metadata with only blank samples;
        if (options.skip <= 4) {
            shell.rm('-rf', __dirname+'/test/L754_bacs/L754_bacs_blanks');
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 4 - L754_bacs_blanks - run @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' run -n L754_bacs_blanks ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata_blanks.csv ' +
                '-d '+__dirname+'/test/L754_bacs/mzxml ' +
                '-o '+__dirname+'/test/L754_bacs/ -y processed_data_blanks -j '+options.tremolo+' -v 11 -q '+
                options.pre_process + ' --bflag_cutoff 1',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[3] = test_res[3] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[3] = test_res[3] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }
            //console.log("\n\n");
        }

        // # test metadata with only one samples;
        if (options.pre_process == "TRUE" && options.skip <= 5) {
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 5 - L754_bacs_blanks_one_sample - pre_process @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec=shell.exec(np3_js_call+' pre_process -n L754_bacs_blanks_one_sample ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata_one_sample.csv ' +
                '-d '+__dirname+'/test/L754_bacs/mzxml -y processed_data_blanks_one_sample -v 11 -q '+options.pre_process,
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[4] = test_res[4] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[4] = test_res[4] + '\n\nDone! :)';
                console.log('DONE!\n');
            }
            //console.log("\n\n");
        } else if (options.skip <= 5) {
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 5 - L754_bacs_blanks_one_sample - pre_process @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\nskipped (pre-processing rerun is disabled)\n");

            test_res[4] = test_res[4] + '\n\nskipped :) (pre-processing rerun is disabled)';
        }

        if (options.skip <= 6) {
            shell.rm('-rf', __dirname+'/test/L754_bacs/L754_bacs_blanks_one_sample');
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 6 - L754_bacs_blanks_one_sample - clustering @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' clustering -n L754_bacs_blanks_one_sample ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata_one_sample.csv ' +
                '-d '+__dirname+'/test/L754_bacs/mzxml ' +
                '-o '+__dirname+'/test/L754_bacs -y processed_data_blanks_one_sample -v 13 -b 500 ' +
                '-q '+options.pre_process+
                ' -j '+options.tremolo,
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[5] = test_res[5] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[5] = test_res[5] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }
            //console.log("\n\n");
        }

        if (options.skip <= 7 && !isWindows()) {
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 7 - L754_bacs_blanks_one_sample - tremolo @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' tremolo ' +
                '-o '+__dirname+'/test/L754_bacs/L754_bacs_blanks_one_sample/outs/L754_bacs_blanks_one_sample/identifications/ ' +
                '-g '+__dirname+'/test/L754_bacs/L754_bacs_blanks_one_sample/outs/L754_bacs_blanks_one_sample/mgf/L754_bacs_blanks_one_sample_all.mgf ' +
                '-k 20 -v 13',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.toLocaleUpperCase().includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[6] = test_res[6] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[6] = test_res[6] + '\n\nDONE :)';
                console.log('DONE!\n');
            }
        } else if (options.skip <= 7) {
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 7 - L754_bacs_blanks_one_sample - tremolo @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\nskipped (tremolo is only available for unix OS)\n");

            test_res[6] = test_res[6] + '\n\nskipped :) (tremolo is only available for unix OS)';
        }

        if (options.skip <= 8) {
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 8 - L754_bacs_blanks_one_sample - clean & annotate_protonated @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' clean ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata_one_sample.csv ' +
                '-o '+__dirname+'/test/L754_bacs/L754_bacs_blanks_one_sample/outs/L754_bacs_blanks_one_sample ' +
                '-y '+__dirname+'/test/L754_bacs/mzxml/processed_data_blanks_one_sample -b 500 -v 13 ' +
                '-t 2,5 --bflag_cutoff 3 --max_shift 500 --min_matched_peaks 10',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.toLocaleUpperCase().includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[7] = test_res[7] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[7] = test_res[7] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }

            resExec = shell.exec(np3_js_call+' annotate_protonated ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata_one_sample.csv ' +
                '-o '+__dirname+'/test/L754_bacs/L754_bacs_blanks_one_sample/outs/L754_bacs_blanks_one_sample ' +
                '-b 500 -v 13 -t 2 -i 500',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.toLocaleUpperCase().includes("ERROR") || resExec.stderr.includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[7] = test_res[7] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[7] = test_res[7] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }
        }

        if (options.skip <= 9) {
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 9 - L754_bacs_blanks_one_sample - merge all @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' merge ' +
                '-o '+__dirname+'/test/L754_bacs/L754_bacs_blanks_one_sample/outs/L754_bacs_blanks_one_sample ' +
                '-y '+__dirname+'/test/L754_bacs/mzxml/processed_data_blanks_one_sample ' +
                '-m '+__dirname+'/test/L754_bacs/marine_bacteria_library_L754_metadata_one_sample.csv -v 13 -p FALSE',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.toLocaleUpperCase().includes("ERROR") || resExec.stderr.toLocaleUpperCase().includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[8] = test_res[8] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[8] = test_res[8] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }
        }

        if (options.skip <= 10) {
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            console.log("@@@@@ Test 10 - L754_bacs_blanks_one_sample - mn @@@@@");
            console.log("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            resExec = shell.exec(np3_js_call+' mn ' +
                '-o '+__dirname+'/test/L754_bacs/L754_bacs_blanks_one_sample/outs/L754_bacs_blanks_one_sample ' +
                '-w 0.9 -k 5 -b 10 -v 13 --min_matched_peaks 1',
                {async:false, silent:true});
            if (resExec.code || resExec.stdout.toLocaleUpperCase().includes("ERROR") || resExec.stderr.toLocaleUpperCase().includes("ERROR")) {
                // in case of error show all the emmited msgs
                console.log(resExec.stdout);
                console.log(resExec.stderr);
                test_res[9] = test_res[9] + '\n\nEXEC ERROR';
                console.log('ERROR\n');
            } else {
                test_res[9] = test_res[9] + resExec.stdout.split('*** TESTING ***\n\n')[1];
                console.log('DONE!\n');
            }
            console.log("\n\n");
        }

        console.log("-------------------------------------------------------\n");
        console.log("--------------------TEST RESULTS-----------------------\n");
        console.log("-------------------------------------------------------\n\n");

        unit_test_res.forEach(function (res) {
            console.log(res+ "\n");
        });
        test_res.forEach(function (res) {
            console.log(res+ "\n");
        });

        console.log("\n-------------------------------------------------------\n");
        console.log(printTimeElapsed_bigint(start_test, process.hrtime.bigint()));
        console.log("\n--------------------------END--------------------------\n");
    });

// error on unknown commands
program.on('command:*', function () {
    console.error('Invalid command: %s\nSee --help for a list of available commands.', program.args.join(' '));
    process.exit(1);
});

program.parse(process.argv);
