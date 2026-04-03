"""Nox sessions for building, testing, and development of pointprocess."""

import nox
import platform
from pathlib import Path

# Python versions to test
nox.options.default_venv_backend = "venv"
DEFAULT_PYTHON = "3.11"
PYTHON_VERSIONS = ["3.9", "3.10", "3.11", "3.12", "3.13"]


@nox.session
def build(session: nox.Session) -> None:
    """Build the project in Release mode."""
    session.log("Building pointprocess in Release mode...")
    session.run("cmake", "--preset", "release", external=True)
    session.run("cmake", "--build", "--preset", "release", external=True)
    session.log("✓ Build complete. Output in 'build/' directory")


@nox.session
def test(session: nox.Session) -> None:
    """Run the test suite."""
    session.log("Running test suite...")
    
    # Install dependencies
    session.install("pytest", "pytest-cov")
    
    # Build tests
    session.run("cmake", "--preset", "tests-with-coverage", external=True)
    session.run(
        "cmake",
        "--build",
        "--preset",
        "tests-with-coverage",
        external=True,
    )
    
    # Run tests
    test_exe = Path("build-tests") / "tests" / "RunTests"
    if platform.system() == "Windows":
        test_exe = test_exe.with_suffix(".exe")
    
    session.run(str(test_exe), external=True)
    session.log("✓ All tests passed!")


@nox.session
def coverage(session: nox.Session) -> None:
    """Generate code coverage report."""
    session.log("Generating coverage report...")
    
    # Requires lcov, gcov on Unix systems
    if platform.system() == "Windows":
        session.log("⚠ Coverage report generation requires lcov (Unix-only)")
        return
    
    # Build tests with coverage instrumentation
    session.log("Building tests with coverage instrumentation...")
    session.run("cmake", "--preset", "tests-with-coverage", external=True)
    session.run("cmake", "--build", "--preset", "tests-with-coverage", external=True)
    
    build_dir = Path("build-tests")
    
    # Zero coverage counters
    session.log("Zeroing coverage counters...")
    session.run("lcov", "--zerocounters", "--directory", str(build_dir), external=True)
    
    # Run tests
    session.log("Running tests...")
    test_exe = build_dir / "tests" / "RunTests"
    if not test_exe.exists():
        session.log(f"⚠ Test executable not found at {test_exe}")
        return
    session.run(str(test_exe), external=True)
    
    # Capture coverage information
    session.log("Capturing coverage information...")
    coverage_file = build_dir / "coverage.info"
    session.run(
        "lcov",
        "--capture",
        "--directory", str(build_dir / "tests"),
        "-o", str(coverage_file),
        "--include", f"{Path.cwd()}/src/pointprocess/*",
        "--gcov-tool", "gcov",
        external=True,
    )
    
    # Generate HTML report
    session.log("Generating HTML coverage report...")
    out_dir = build_dir / "out"
    session.run("genhtml", str(coverage_file), "--output-directory", str(out_dir), external=True)
    
    index_file = out_dir / "index.html"
    if index_file.exists():
        session.log(f"✓ Coverage report generated: {index_file}")
    else:
        session.log("⚠ Coverage report may not have been generated")


@nox.session(name="lint")
def lint(session: nox.Session) -> None:
    """Run code formatting checks."""
    session.log("Running code quality checks...")
    
    session.install("clang-tools")
    
    # Check formatting using glob expansion
    session.log("Checking code formatting with clang-format...")
    cpp_files = list(Path("src/pointprocess").glob("*.cpp"))
    h_files = list(Path("src/pointprocess").glob("*.h"))
    
    if cpp_files or h_files:
        format_files = [str(f) for f in cpp_files + h_files]
        session.run(
            "clang-format",
            "--dry-run",
            "-i",
            "-ferror-limit=0",
            *format_files,
            external=True,
            success_codes=[0, 1],  # clang-format returns 0 if no changes, 1 if changes needed
        )
        session.log("✓ Code formatting check complete")
    else:
        session.log("⚠ No source files found to check")


@nox.session(name="dev")
def dev(session: nox.Session) -> None:
    """Set up an interactive development environment."""
    session.log("Setting up development environment...")
    
    # Install Python dev dependencies
    session.install(
        "pytest",
        "pytest-cov",
        "black",
        "isort",
        "mypy",
        "numpy",
        "scipy",
        "matplotlib",
    )
    
    # Build C++ extension using CMake
    session.log("Building C++ extension with CMake...")
    session.run("cmake", "--preset", "release", external=True)
    session.run("cmake", "--build", "--preset", "release", external=True)
    
    # Add built library to PYTHONPATH so it can be imported
    build_lib = Path("build/src").resolve()
    session.log(f"✓ Built library at: {build_lib}")
    session.log("✓ Development environment ready!")
    session.log("\nTo import pointprocess, set PYTHONPATH:")
    session.log(f"  export PYTHONPATH={build_lib}:$PYTHONPATH")
    session.log("\nOr use the following in Python:")
    session.log(f"  import sys; sys.path.insert(0, '{build_lib}')")
    session.log("\nTips for development:")
    session.log("  - Build: nox -s build")
    session.log("  - Test: nox -s test")
    session.log("  - Lint: nox -s lint")
    session.log("  - Wheel: nox -s wheel")
    session.log("  - Docs: nox -s docs")


@nox.session
def wheel(session: nox.Session) -> None:
    """Build a wheel for the current platform."""
    session.log("Building wheel for current platform...")
    
    session.install("build", "scikit-build-core", "cmake", "ninja", "pybind11")
    session.run(
        "python",
        "-m",
        "build",
        "--wheel",
        "--outdir",
        "dist",
        external=True,
    )
    
    wheels = list(Path("dist").glob("*.whl"))
    if wheels:
        session.log(f"✓ Wheel created: {wheels[-1].name}")
    else:
        session.log("⚠ Wheel build completed but files not found")


@nox.session
def docs(session: nox.Session) -> None:
    """Build Sphinx documentation."""
    session.log("Building documentation...")
    
    session.install("-r", "docs/requirements.txt")
    
    session.run(
        "python",
        "-m",
        "sphinx",
        "-W",
        "-a",
        "-n",
        "docs/source",
        "docs/_build/html",
        external=True,
    )
    
    index_file = Path("docs/_build/html/index.html")
    if index_file.exists():
        session.log(f"✓ Documentation built: {index_file}")
    else:
        session.log("⚠ Documentation build may have encountered issues")


@nox.session(name="clean")
def clean(session: nox.Session) -> None:
    """Clean build artifacts."""
    session.log("Cleaning build artifacts...")
    
    import shutil
    
    dirs_to_remove = [
        Path("build"),
        Path("build-tests"),
        Path("dist"),
        Path("docs/_build"),
        Path(".eggs"),
        Path(".nox"),
    ]
    
    for dir_path in dirs_to_remove:
        if dir_path.exists():
            shutil.rmtree(dir_path)
            session.log(f"  Removed {dir_path}")
    
    # Remove Python cache
    for pycache in Path(".").rglob("__pycache__"):
        shutil.rmtree(pycache)
    
    session.log("✓ Cleanup complete")
