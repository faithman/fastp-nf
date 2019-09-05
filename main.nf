#!/usr/bin/env nextflow

params.email = 'yewangfaith@gmail.com'
params.directory = "$PWD/"
params.out = params.directory.replace("raw", "test")
println params.out
params.threads = 8
println "Running Fastp trim on " + params.directory
println params.directory + '*_{1,2}.fq.gz'

// Fetch fqs; alternative suffixes
Channel.fromFilePairs(params.directory + '*_{1,2}.fq.gz', flat: true)
        .into { pre_trim_fastqc; trimmomatic_read_pairs; pre_trim_fastqc; md5_fastqc; log_fq }
        
log_fq.subscribe { println it }

process make_out_dir {
    
    executor 'local'
    
    """
    mkdir -p ${params.out}
    """
}

/* 
    ======================
    Perform initial FASTQC
    ======================
*/

process pre_trim_fastqc {

    tag { dataset_id } 
    
    publishDir "${params.directory}/fastqc", mode: 'copy'
        
    cpus 8
    
    input:
        set dataset_id, file(forward), file(reverse) from pre_trim_fastqc
    
    output:
        set dataset_id, file("${dataset_id}_1_fastqc.zip"), file("${dataset_id}_1_fastqc.zip"), file("${dataset_id}_2_fastqc.html"), file("${dataset_id}_2_fastqc.html") into pre_trim_multi_qc
    """
        fastqc --noextract --threads ${task.cpus} ${forward}
        fastqc --noextract --threads ${task.cpus} ${reverse}
    """
}

/* 
    ============================
    Aggregate pre-FASTQC results
    ============================
*/

process pre_trim_multi_qc_run {

    publishDir "${params.out}/report", mode: 'copy'

    input:
        file(dataset_id) from pre_trim_multi_qc.toSortedList()
    output:
        file("multiqc_report_pre.html")
    """
        multiqc --filename multiqc_report_pre.html ${params.directory}/fastqc
    """

}

/* 
    =======================
    Generate an MD5sum hash
    =======================
*/

process md5sum_pre {

    tag { dataset_id } 

    input:
        set dataset_id, file(forward), file(reverse) from md5_fastqc
    output:
        file('md5.txt') into md5_set
    """
    # Command to determine md5
    
    __rvm_md5_for()
    {
      if builtin command -v md5 > /dev/null; then
        echo \"\$1\" | md5
      elif builtin command -v md5sum > /dev/null ; then
        echo \"\$1\" | md5sum | awk '{print \$1}'
      else
        rvm_error "Neither md5 nor md5sum were found in the PATH"
        return 1
      fi

      return 0
    }

        __rvm_md5_for ${forward} | awk '{ print \$0 "\\t${dataset_id}\\t${forward}" }' > md5.txt
        __rvm_md5_for ${reverse} | awk '{ print \$0 "\\t${dataset_id}\\t${reverse}" }' >> md5.txt
    """
}

process collect_md5sum {

    publishDir "${params.out}/report", mode: 'copy'

    input:
        file(md5) from md5_set.collectFile(name: 'md5sum.txt')

    output:
        file("md5sum.txt")

    """
        echo "Great!"
    """

}

/* 
    ================
    Perform trimming
    ================
*/

process trim {

    publishDir params.out, mode: 'copy'

    cpus 8
    
    tag { dataset_id } 
    
    input:
        set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

    output:
        set dataset_id, file("${dataset_id}_1P.fq.gz"), file("${dataset_id}_2P.fq.gz") into trim_output

    """
    fastp -t 4 --detect_adapter_for_pe -i $forward -I $reverse -o ${dataset_id}_1P.fq.gz -O ${dataset_id}_2P.fq.gz
    """

}

/* 
    ===========
    Post FASTQC
    ===========
*/


process post_trim_fastqc {

    publishDir "${params.out}/fastqc", mode: 'copy'
    
    stageInMode 'copy'
    
    cpus 8
    
    tag { dataset_id } 
    
    input:
        set dataset_id, file("${dataset_id}_1P.fq.gz"), file("${dataset_id}_2P.fq.gz") from trim_output
    
    output:
        set file("${dataset_id}_1P_fastqc.zip"), file("${dataset_id}_2P_fastqc.zip"), file("${dataset_id}_1P_fastqc.html"), file("${dataset_id}_2P_fastqc.html") into post_trim_multi_qc
    
    """
        fastqc --noextract --threads ${task.cpus} ${dataset_id}_1P.fq.gz
        fastqc --noextract --threads ${task.cpus} ${dataset_id}_2P.fq.gz
    """
}

process post_trim_multi_qc_run {

    publishDir "${params.out}/report", mode: 'copy'

    input:
        file(dataset_id) from post_trim_multi_qc.toSortedList()
    
    output:
        file("multiqc_report_post.html")
        
    """
        multiqc --filename multiqc_report_post.html ${params.out}/fastqc
    """

}


workflow.onComplete {

    user="whoami".execute().text

    summary = """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    User: ${user}
    """

    println summary

    // mail summary
    ['mail', '-s', 'wi-nf', params.email].execute() << summary

    def outlog = new File("${params.out}/report/trimming_log.txt")
    outlog.newWriter().withWriter {
        outlog << summary
        outlog << "\n--------pyenv-------\n"
        outlog << "pyenv versions".execute().text
        outlog << "--------ENV--------"
        outlog << "ENV".execute().text
        outlog << "--------brew--------"
        outlog << "brew list".execute().text
        outlog << "--------R--------"
        outlog << "Rscript -e 'devtools::session_info()'".execute().text
    }

}

