"""
Implementation of the CIME LILACSMOKE (LILAC smoke) test.

This is a CTSM-specific test. It tests the building and running of CTSM via LILAC. Compset
is ignored, but grid is important. Also, it's important that this test use the nuopc
driver, both for the sake of the build and for extracting some runtime settings. This test
should also use the lilac testmod (or a testmod that derives from it) in order to
establish the user_nl_ctsm file correctly.
"""

import os
import shutil

from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.utils import run_cmd_no_fail, append_testlog, symlink_force
from CIME.build import post_build
from CIME.test_status import GENERATE_PHASE, BASELINE_PHASE, TEST_PASS_STATUS
from CIME.XML.standard_module_setup import *

logger = logging.getLogger(__name__)

class LILACSMOKE(SystemTestsCommon):

    def __init__(self, case):
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):
        if not sharedlib_only:
            caseroot = self._case.get_value('CASEROOT')
            lndroot = self._case.get_value('COMP_ROOT_DIR_LND')
            exeroot = self._case.get_value('EXEROOT')
            build_dir = os.path.join(caseroot, 'lilac_build')
            script_path = os.path.abspath(os.path.join(lndroot, 'lilac', 'build_ctsm'))

            # We only run the initial build command if the build_dir doesn't exist
            # yet. This is to support rebuilding the test case. (The first time through,
            # the build_dir won't exist yet; subsequent times, it will already exist, so
            # we skip to the rebuild command.)
            if not os.path.isdir(build_dir):
                machine = self._case.get_value('MACH')
                compiler = self._case.get_value('COMPILER')
                debug = self._case.get_value('DEBUG')
                # It would be possible to do this testing via the python interface rather
                # than through a separate subprocess. However, we do it through a
                # subprocess in order to test the full build_ctsm script, including
                # command-line parsing.
                cmd = '{script_path} {build_dir} --machine {machine} --compiler {compiler}'.format(
                    script_path=script_path,
                    build_dir=build_dir,
                    machine=machine,
                    compiler=compiler)
                if debug:
                    cmd += ' --build-debug'
                self._run_build_cmd(cmd, exeroot, 'build_ctsm.bldlog')

            # We call the build script with --rebuild even for an initial build. This is
            # so we make sure to test the code path for --rebuild. (This is also needed if
            # the user rebuilds the test case, in which case this will be the only command
            # run, since the build_dir will already exist.)
            cmd = '{script_path} {build_dir} --rebuild'.format(
                script_path=script_path,
                build_dir=build_dir)
            self._run_build_cmd(cmd, exeroot, 'rebuild_ctsm.bldlog')

            self._build_atm_driver()

            self._create_runtime_inputs()

            self._setup_atm_driver_rundir()

            # Setting logs=[] implies that we don't bother gzipping any of the build log
            # files; that seems fine for these purposes (and it keeps the above code
            # simpler).
            post_build(self._case, logs=[], build_complete=True)

    def _build_atm_driver(self):
        caseroot = self._case.get_value('CASEROOT')
        lndroot = self._case.get_value('COMP_ROOT_DIR_LND')
        blddir = os.path.join(caseroot, 'lilac_atm_driver', 'bld')

        if not os.path.exists(blddir):
            os.makedirs(blddir)
        symlink_force(os.path.join(lndroot, 'lilac', 'atm_driver', 'Makefile'),
                      os.path.join(blddir, 'Makefile'))
        symlink_force(os.path.join(lndroot, 'lilac', 'atm_driver', 'atm_driver.F90'),
                      os.path.join(blddir, 'atm_driver.F90'))
        symlink_force(os.path.join(caseroot, 'Macros.make'),
                      os.path.join(blddir, 'Macros.make'))

        makevars = 'COMPILER={compiler} DEBUG={debug} CTSM_MKFILE={ctsm_mkfile}'.format(
            compiler=self._case.get_value('COMPILER'),
            debug=str(self._case.get_value('DEBUG')).upper(),
            ctsm_mkfile=os.path.join(caseroot, 'lilac_build', 'ctsm.mk'))
        makecmd = 'make {makevars} atm_driver'.format(makevars=makevars)
        self._run_build_cmd(makecmd, blddir, 'atm_driver.bldlog')

    def _create_runtime_inputs(self):
        caseroot = self._case.get_value('CASEROOT')
        runtime_inputs = self._runtime_inputs_dir()
        lnd_domain_file = os.path.join(self._case.get_value('LND_DOMAIN_PATH'),
                                       self._case.get_value('LND_DOMAIN_FILE'))

        # Cheat a bit here: Get the fsurdat file from the already-generated lnd_in file in
        # the host test case - i.e., from the standard cime-based preview_namelists. But
        # this isn't really a morally-objectionable cheat, because in the real workflow,
        # we expect the user to identify fsurdat manually; in this testing situation, we
        # need to come up with some way to replace this manual identification, so cheating
        # feels acceptable.
        self._case.create_namelists(component='lnd')
        fsurdat = self._extract_var_from_namelist(
            nl_filename=os.path.join(caseroot, 'CaseDocs', 'lnd_in'),
            varname='fsurdat')

        self._fill_in_variables_in_file(filepath=os.path.join(runtime_inputs, 'ctsm.cfg'),
                                        replacements={'lnd_domain_file':lnd_domain_file,
                                                      'fsurdat':fsurdat})

        # The user_nl_ctsm in the case directory is set up based on the standard testmods
        # mechanism. We use that one in place of the standard user_nl_ctsm, since the one
        # in the case directory may contain test-specific modifications.
        shutil.copyfile(src=os.path.join(caseroot, 'user_nl_ctsm'),
                        dst=os.path.join(runtime_inputs, 'user_nl_ctsm'))

        self._run_build_cmd('make_runtime_inputs --rundir {}'.format(runtime_inputs),
                            runtime_inputs,
                            'make_runtime_inputs.log')

        # In lilac_in, we intentionally use the land mesh file for both atm_mesh_filename
        # and lnd_mesh_filename
        lnd_mesh = self._case.get_value('LND_DOMAIN_MESH')
        self._fill_in_variables_in_file(filepath=os.path.join(runtime_inputs, 'lilac_in'),
                                        replacements={'atm_mesh_filename':lnd_mesh,
                                                      'lnd_mesh_filename':lnd_mesh})

        # We run download_input_data partly because it may be needed and partly to test
        # this script. Note, though, that some files (the surface dataset, land domain
        # file, and land/atm mesh file) will not be downloaded automatically: Because of
        # the way this test fills in these file paths, they end up pointing to the true
        # inputdata directory on this machine rather than the sym linked one in
        # lilac_build_dir (which kind of makes sense, because this simulates files that
        # the user would provide themselves). Thus, these files live outside of the
        # inputdata root considered by download_input_data. So if these files are missing,
        # they will need to be downloaded manually (e.g., using check_input_data from the
        # test case).
        self._verify_inputdata_link()
        self._run_build_cmd('download_input_data --rundir {}'.format(runtime_inputs),
                            runtime_inputs,
                            'download_input_data.log')

    def _setup_atm_driver_rundir(self):
        caseroot = self._case.get_value('CASEROOT')
        lndroot = self._case.get_value('COMP_ROOT_DIR_LND')
        rundir = os.path.join(caseroot, 'lilac_atm_driver', 'run')

        if not os.path.exists(rundir):
            os.makedirs(rundir)
            shutil.copyfile(src=os.path.join(lndroot, 'lilac', 'atm_driver', 'atm_driver_in'),
                            dst=os.path.join(rundir, 'atm_driver_in'))

        # As above: assume the land variables also apply to the atmosphere
        lnd_mesh = self._case.get_value('LND_DOMAIN_MESH')
        lnd_nx = self._case.get_value('LND_NX')
        lnd_ny = self._case.get_value('LND_NY')
        expect(self._case.get_value('STOP_OPTION') == 'ndays',
               'LILAC testing currently assumes STOP_OPTION of ndays, not {}'.format(
                   self._case.get_value('STOP_OPTION')))
        stop_n = self._case.get_value('STOP_N')
        self._fill_in_variables_in_file(filepath=os.path.join(rundir, 'atm_driver_in'),
                                        replacements={'atm_mesh_file':lnd_mesh,
                                                      'atm_global_nx':str(lnd_nx),
                                                      'atm_global_ny':str(lnd_ny),
                                                      'atm_stop_day':str(stop_n+1),
                                                      'atm_ndays_all_segs':str(stop_n)})

    def _extract_var_from_namelist(self, nl_filename, varname):
        """Tries to find a variable named varname in the given file; returns its value

        If not found, aborts
        """
        with open(nl_filename) as nl_file:
            for line in nl_file:
                match = re.search(r'^ *{} *= *[\'"]([^\'"]+)'.format(varname), line)
                if match:
                    return match.group(1)
        expect(False, '{} not found in {}'.format(varname, nl_filename))

    def _fill_in_variables_in_file(self, filepath, replacements):
        """For the given file, make the given replacements

        replacements should be a dictionary mapping variable names to their values
        """
        orig_filepath = '{}.orig'.format(filepath)
        if not os.path.exists(orig_filepath):
            shutil.copyfile(src=filepath,
                            dst=orig_filepath)
        os.remove(filepath)

        with open(orig_filepath) as orig_file:
            with open(filepath, 'w') as new_file:
                for orig_line in orig_file:
                    line = orig_line
                    for varname in replacements:
                        # In practice, at most one of the replacements will be done. But
                        # it's simplest just to use a loop over all possible replacements.
                        line = self._fill_in_variable(line=line,
                                                      varname=varname,
                                                      value=replacements[varname])
                    new_file.write(line)

    def _fill_in_variable(self, line, varname, value):
        """Fill in a FILL_THIS_IN variable in a config or namelist file

        Returns the line with FILL_THIS_IN replaced with the given value if this line is
        for varname; otherwise returns line unchanged.
        """
        if re.search(r'^ *{} *='.format(varname), line):
            newline = line.replace('FILL_THIS_IN', value)
        else:
            newline = line
        return newline

    def _verify_inputdata_link(self):
        """Verify that the inputdata link has been set up correctly"""
        din_loc_root = self._case.get_value('DIN_LOC_ROOT')
        inputdata = os.path.join(self._case.get_value('CASEROOT'), 'lilac_build', 'inputdata')
        expect(os.path.realpath(inputdata) == os.path.realpath(din_loc_root),
               'inputdata not set up with the correct link')

    def _runtime_inputs_dir(self):
        return os.path.join(self._case.get_value('CASEROOT'), 'lilac_build', 'runtime_inputs')

    @staticmethod
    def _run_build_cmd(cmd, exeroot, logfile):
        """
        Runs the given build command, with output to the given logfile

        Args:
        cmd: str (command to run)
        exeroot: str (path to exeroot)
        logfile: str (path to logfile)
        """
        append_testlog(cmd)
        run_cmd_no_fail(cmd, arg_stdout=logfile, combine_output=True, from_dir=exeroot)
        with open(os.path.join(exeroot, logfile)) as lf:
            append_testlog(lf.read())

    def run_phase(self):
        # FIXME(wjs, 2020-06-10) Fill this in
        pass

