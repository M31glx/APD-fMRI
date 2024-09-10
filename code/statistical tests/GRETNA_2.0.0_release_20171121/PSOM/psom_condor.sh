#! /bin/sh
# Written by Michael Hanke <michael.hanke@gmail.com>
# To generate a manpage, simply run through help2man.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# play safe
set -u
set -e

# define variable defaults
scriptname="condor_qsub"
version="0.1"
error_file=""
binary_flag=0
cwd_flag=0
hold_jid=""
mail_config=""
mail_user=""
job_name=""
output_file=""
priority=0
rerun_flag=0
shell=${SHELL:-/bin/sh}
taskid_first=""
taskid_last=""
taskid_step=""
qsub_cmd=""
qsub_args=""
queue_name=""
export_env_flag=0
condor_keep_files=0

print_description()
{
cat << EOT
Minimalistic emulation of SGE's qsub for Condor.

EOT
}

print_version()
{
cat << EOT
$scriptname $version

This script is released under the Apache V2.0 License. It was
written by Michael Hanke <michael.hanke@gmail.com>.

EOT
}

print_help()
{
cat << EOT
Usage:  $scriptname [OPTIONS] [<command> [<command_args>]]

The primary purpose of this emulation is to allow SGE-style submission of
dependent jobs without the need to specify the full dependency graph at
submission time. This implementation is neither as efficient as Condor's
DAGMan, nor as functional as SGE's qsub/qalter. It merely serves as a minimal
adaptor to be able to use software original written to interact with SGE in a
Condor pool.

In general $scriptname behaves just like qsub. However, only a fraction of the
original functionality is available. The following list of options only
describes the differences in the behavior of SGE's qsub and this emulation.
Qsub options not listed here are not supported.

Options:

-b <y|n>
  If 'y', command and arguments given on the command line are wrapped into a
  shell script which is then submitted to Condor.

--condor-keep-files
  This is a non-SGE option. If given, it will prevent $scriptname from deleting
  temporary files (generated submit files, sentinel jobs). This is mostly useful
  for debugging.

-cwd
  If given, this option will cause the 'initialdir' value in the Condor submit
  file to be set to the current directory.

-e <filename|path>
  Name of the file to contain the STDERR output of the job. By default this will
  be job_name.ejob_id[.task_id]. If an existing directory is specified, the file
  will be placed inside this directory using the default schema for the
  filename.

-h,--help
  Print usage summary and option list.

-hold_jid <jid>
  If given, the job will be submitted in 'hold' state. Along with the actual job
  a 'sentinel' job will be submitted to Condor's local universe. This sentinel
  watches the specified job and releases the submitted job whenever the job has
  completed. The sentinel observes SGE's behavior to detect job exiting with
  code 100 and not start depedent job in this case. If a cluster id of an array
  job is given the dependent job will only be released after all individual jobs
  of a cluster have completed.

-l <ressource spec>
  This option is currently ignored.

-m <a|e|n><...>
  SGE's notification labels will be translated (approximately) into Condor's
  notifications states (Never, Error, Complete).

-M <email>
  Added as 'notify_user' to the submit file.

-N <jobname>
  Determines the default name of logfile (stdout, stderr).

-o <filename|path>
  See -e option, but for a job's stdout.

-p <int>
  Added a 'priority' to the submit file.

-r <y|n>
  This option is currently ignored.

-S <shell>
  Path to a shell binary for script execution.

-shell <y|n>
  This option is currently ignored.

-t <start>[-<stop>[:<step>]]
  Task ID specification for array job submissions.

-q <queue name>
  This option is permanently ignored, as Condor doesn't have multiple queues.

-V
  If given, 'getenv = True' is added to the submit file.

--version
  Print version information and exit.

EOT
}

create_dep_sentinel() {
	cat << EOT > $1
#!/bin/sh

hold_jids="$2"
dep_job="$3"

clean_error() {
	printf "\$1\n"
	condor_rm \$dep_job
	exit 1
}

# as long as there are relevant job in the queue wait and try again
while [ \$(condor_q -long -attributes Owner \$hold_jids | grep "$USER" | wc -l) -ge 1 ]; do
	sleep 5
done

# now check whether all deps have exited properly (i.e. not with code 100)
for jid in \$hold_jids; do
	   job_status="\$(condor_history \$jid |tail -n 1 | awk '{print \$6}')"
	   case "\$job_status" in
		   C) if [ "\$(condor_history -long \$jid |grep ExitCode | cut -d= -f2,2)" = " 100" ]; \
			  then clean_error "Error: Job dependency \$jid completed but exited with code 100."; fi;;
			X) clean_error "Error: Job dependency \$jid has been removed from the queue" ;;
			*) clean_error "Error: Job dependency \$jid doesn't exist" ;;
		esac
done

# all good -- let the job go
condor_release \$dep_job

# in the end we can safely clean this sentinel script
[ $condor_keep_files -eq 0 ] && rm \$0 || true

EOT
	chmod +x $1

cat << EOT > $1.submit
universe = local
executable = $1
Queue
EOT
	chmod +x $1.submit

}

parse_args() {
	###############################################################
	# cmdline args
	###############################################################

	# Parse commandline options
	# Note that we use `"$@"' to let each command-line parameter expand to a
	# separate word. The quotes around `$@' are essential!
	# We need CLOPTS as the `eval set --' would nuke the return value of getopt.
	CLOPTS=`getopt -a -o b:,e:,h,l:,m:,M:,N:,o:,p:,r:,S:,t:,q:,V --long cwd,help,hold_jid:,shell:,verbose,version,condor-keep-files -n "$scriptname" -- "$@"`

	if [ $? != 0 ] ; then
	  echo "Terminating..." >&2
	  exit 1
	fi

	# Note the quotes around `$CLOPTS': they are essential!
	eval set -- "$CLOPTS"

	while true ; do
	  case "$1" in
		-b) shift; if [ "$1" = "y" ]; then binary_flag=1; fi; shift;;
		--cwd) shift; cwd_flag=1;;
		--condor-keep-files) shift; condor_keep_files=1;;
		-e) shift; error_file=$1; shift;;
		-h|--help) print_description; print_help; exit 0;;
		--hold_jid) shift; hold_jid=$1; shift;;
		# we would be only interested in the arch spec -- right now it is ignored
		-l) shift; shift;;
		-m) shift; mail_config=$1; shift;;
		-M) shift; mail_user=$1; shift;;
		-N) shift; job_name=$1; shift;;
		-o) shift; output_file=$1; shift;;
		-p) shift; priority=$1; shift;;
		-r) shift; rerun_flag=$1; shift;;
		-S) shift; shell=$1; shift;;
		--shell) shift; shift;;
		-t) shift
			taskid_first="$(echo "$1" | awk -F- '{print $1}')"
			taskid_last="$(echo "$1" | awk -F- '{print $2}' | awk -F: '{print $1}')"
			taskid_step="$(echo "$1" | awk -F- '{print $2}' | awk -F: '{print $2}')"
			shift;;
		# needs to handle SGE_TASK_ID, SGE_TASK_FIRST, SGE_TASK_LAST and SGE_TASK_STEPSIZE
		# log goes to: <jobname>.['e'|'o']<job_id>'.'<task_id>
		-q) shift; queue_name=$1; shift;;
		-V) shift; export_env_flag=1;;
		--version) print_version; exit 0;;
		--) shift ; break ;;
		*) echo "Internal error! ($1)"; exit 1;;
	  esac
	done
	# the first arg is the command the rest are its arguments
	if [ $# -ge 1 ] && [ "$1" != "XXcondor_sub_scriptmodeXX" ]; then
		qsub_cmd="$1"
		shift
		qsub_args="$@"
	fi
}

# parse all commandline args
parse_args $@

# no arguments -> need to read from stdin
if [ -z "$qsub_cmd" ]; then
	# redirect STDIN into a file, place it into the current dir to increase the
	# likelihood of being available on the execute host too (NFS-mounted FS)
	# unfortunately the author cannot think of a way to clean this tempfile up
	# as it need to be available until the last job in the cluster actually
	# started running -- yet another sentinel could do it ....
	if [ -z "$job_name" ]; then
		stdin_script_file=$(mktemp --tmpdir=$(pwd) STDIN_qsub.XXXXXXXXXXXXX)
	else
		stdin_script_file=$(mktemp --tmpdir=$(pwd) ${job_name}_qsub.XXXXXXXXXXXXX)
	fi
	cat < /dev/stdin > $stdin_script_file
	chmod +x $stdin_script_file
	# use same default name as SGE
	if [ -z "$job_name" ]; then job_name="STDIN"; fi
	qsub_cmd="$stdin_script_file"
	qsub_args=""
fi

# if we're not in "binary" mode, we also need to parse the submitted script for
# additional arguments -- in an ideal world this would be done before parsing
# the actual commandline args to give them a higher priority
if [ $binary_flag -ne 1 ]; then
	script_args=`grep '^\#\$ ' < $qsub_cmd |cut -d' ' -f 2- | tr "\n" " "`
	if [ $? = 0 ] ; then
		# found some args
		parse_args $script_args XXcondor_sub_scriptmodeXX
	fi
fi

# have a meaningful job name in any case
if [ -z "$job_name" ]; then
	job_name="$(basename $qsub_cmd)"
fi

# handle job dependency
# store effective job deps
job_deps=""
if [ -n "$hold_jid" ]; then
	# loop over potentially multiple job deps
	for jid in $(echo "$hold_jid" |tr "," " ") ; do
		if [ $(condor_q -long -attributes Owner $jid | grep "$USER" | wc -l) -ge 1 ]; then
			# job is currently in the queue and owned by $USER -> store
			job_deps="$job_deps $jid"
		else
			# job not owned by this user or not in the queue, but maybe already
			# done?
		   job_status="$(condor_history $jid |tail -n 1 | awk '{print $6}')"
		   echo "JOBSTATUS: $job_status" >&2
		   case "$job_status" in
			   C) if [ "$(condor_history -long $jid |grep ExitCode | cut -d= -f2,2)" = " 100" ]; \
			      then printf "Error: Job dependency $jid completed but exited with code 100.\n"; exit 100; fi;;
				X) printf "Error: Job dependency $jid has been removed from the queue\n" ; exit 1;;
				*) printf "Error: Job dependency $jid doesn't exist\n" ; exit 1;;
			esac
		fi
	done
fi

# compose the names of logfiles
log_file="${log_file:-$job_name.log\$(Cluster)}"

# SGE also accepts directories for stdout and stderr. in this case it would
# place a file with the default name in this directory. This is probably
# evaluated at execution time. for condor we can only easily do this at
# submission time -- this is not the same -- let's hope it works nevertheless

if [ -d "$error_file" ]; then
	error_file="${error_file}/$job_name.e\$(Cluster)"
else
	error_file="${error_file:-$job_name.e\$(Cluster)}"
fi
if [ -d "$output_file" ]; then
	output_file="${output_file}/$job_name.o\$(Cluster)"
else
	output_file="${output_file:-$job_name.o\$(Cluster)}"
fi

# main submit file
submit_file=$(mktemp --tmpdir condor_qsub.XXXXXXXXXXXXX)

# write submit file header
cat << EOT > $submit_file
# condor_qsub call: $@
universe = vanilla
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
#log = $log_file
EOT

# handle 'binary mode'
if [ $binary_flag = 1 ]; then
	## go with current shell or POSIX if undefined
	cat << EOT >> $submit_file
executable = $shell
arguments = "-c '$qsub_cmd $qsub_args'"
EOT
	qsub_cmd="$shell"
else
	cat << EOT >> $submit_file
executable = $shell
arguments = $qsub_cmd $qsub_args
EOT
fi

if [ $cwd_flag = 1 ]; then
	printf "initialdir = $(pwd)\n" >> $submit_file
fi
if [ -n "$error_file" ]; then
	printf "error = $error_file\n" >> $submit_file
fi
if [ -n "$output_file" ]; then
	printf "output = $output_file\n" >> $submit_file
fi
if [ $export_env_flag = 1 ]; then
	printf "getenv = True\n" >> $submit_file
fi
if [ -n "$mail_user" ]; then
	printf "notify_user = $mail_user\n" >> $submit_file
fi
if [ -n "$mail_config" ]; then
	if [ "$mail_config" != "$(echo "$mail_config" | sed -e 's/n//')" ]; then
		printf "notification = Never\n" >> $submit_file
	elif [ "$mail_config" != "$(echo "$mail_config" | sed -e 's/e//')" ]; then
		printf "notification = Complete\n" >> $submit_file
	elif [ "$mail_config" != "$(echo "$mail_config" | sed -e 's/a//')" ]; then
		printf "notification = Error\n" >> $submit_file
	fi
else
	printf "notification = Never\n" >> $submit_file
fi
if [ -n "$priority" ]; then
	printf "priority = $priority\n" >> $submit_file
fi
if [ -n "$job_deps" ]; then
	# in case of job deps, submit the job held and release it via a sentinel
	printf "hold = True\n" >> $submit_file
fi

# handle array jobs
if [ -n "$taskid_first" ]; then
	for taskid in $(seq $taskid_first $taskid_step $taskid_last); do
		printf "environment = \"SGE_TASK_ID=$taskid SGE_TASK_FIRST=$taskid_first SGE_TASK_STEPSIZE=$taskid_step SGE_TASK_LAST=$taskid_last\"\n" >> $submit_file
		if [ -n "$error_file" ]; then
			printf "error = $error_file.$taskid\n" >> $submit_file
		fi
		if [ -n "$output_file" ]; then
			printf "output = $output_file.$taskid\n" >> $submit_file
		fi
		# queue this job
		printf "Queue\n" >> $submit_file
	done
else
	# unconditional queuing
	printf "Queue\n" >> $submit_file
fi

# Done with creating to submission file

# actually submit the job, printing only the cluster id
cluster_id="$(condor_submit - < $submit_file |grep '^** Proc' || true)"
if [ -z "$cluster_id" ]; then
	printf "Job submission failed.\n"
	[ $condor_keep_files -eq 0 ] && rm $submit_file || true
	exit 1
fi
cluster_id="$(echo "$cluster_id" | head -n1 | cut -d ' ' -f 3,3 | cut -d. -f1,1)"

#cleanup
[ $condor_keep_files -eq 0 ] && rm $submit_file || true

# deal with local job sentinel
if [ -n "$job_deps" ]; then
	sentinel_file="$(mktemp --tmpdir cluster${cluster_id}_sentinel.XXXXXXXXXXXXX)"
	create_dep_sentinel $sentinel_file "$job_deps" "$cluster_id"
	condor_submit - < $sentinel_file.submit > /dev/null
	[ $condor_keep_files -eq 0 ] && rm $sentinel_file.submit || true
fi

printf "Your job $cluster_id (\"$job_name\") has been submitted\n"
exit 0

