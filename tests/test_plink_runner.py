from deepg2p_celiac.qc.plink_runner import build_plink_qc_command, run_plink_qc


def test_build_plink_qc_command() -> None:
    command = build_plink_qc_command("input", "output")

    assert command.as_list()[:3] == ["plink", "--bfile", "input"]
    assert "--make-bed" in command.as_list()
    assert command.as_list()[-2:] == ["--out", "output"]


def test_run_plink_qc_dry_run() -> None:
    command = build_plink_qc_command("input", "output")

    assert run_plink_qc(command, dry_run=True) is None

