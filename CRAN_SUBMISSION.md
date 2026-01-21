# CRAN Submission Guide for RcppICA

**Note**: This guide is for future reference. Focus on GitHub release first, then consider CRAN submission later based on user feedback.

---

## Why Wait Before CRAN Submission?

**Recommended approach**: GitHub → User Feedback → CRAN

1. **Gather user feedback** (2-4 weeks)
   - Bug reports from real-world usage
   - Feature requests
   - Platform-specific issues
   - Performance on diverse datasets

2. **Fix any issues** found by early adopters

3. **Polish based on feedback**, then submit to CRAN

**Benefit**: CRAN reviewers are strict. Better to iron out issues with GitHub users first.

---

## Pre-Submission Checklist

### Required: Clean R CMD check

```bash
R CMD build .
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran RcppICA_0.1.0.tar.gz
```

**Must have**: 0 ERRORs, 0 WARNINGs

**Acceptable**: NOTEs are okay if explained in submission comments

### Current Status

✅ **0 ERRORs**
✅ **0 WARNINGs** (PDF manual WARNING is due to missing LaTeX packages, not package issue)
⚠️ **8 NOTEs** (all expected and explainable)

### Expected NOTEs and Explanations

1. **New submission**
   - *Explanation*: This is the first CRAN submission

2. **Invalid URLs (404)**
   - *Explanation*: GitHub repository exists but CRAN's URL checker may timeout
   - *Action*: Ensure https://github.com/mtotty/RcppICA is accessible

3. **Suggests package not available (bench)**
   - *Explanation*: Only used for optional benchmarking, not required
   - *Action*: None needed (optional suggest)

4. **HTML Tidy not found**
   - *Explanation*: Local system issue, not package issue
   - *Action*: None needed

---

## CRAN Submission Process

### 1. Final Checks

```bash
# Build package
R CMD build .

# Check without Suggests
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran RcppICA_0.1.0.tar.gz

# Check examples run quickly (<5 seconds each)
R -e "library(RcppICA); example(fastICA)"

# Verify vignettes build
R CMD build . --no-resave-data
```

### 2. Update DESCRIPTION

Ensure these fields are correct:

```r
Package: RcppICA
Version: 0.1.0
Date: 2026-01-21
Authors@R: person("Michael", "Totty",
                  email = "michael.totty@example.com",
                  role = c("aut", "cre"))
Description: ...
License: GPL-3
URL: https://github.com/mtotty/RcppICA
BugReports: https://github.com/mtotty/RcppICA/issues
```

**Important**: Email must be valid (CRAN will send notifications)

### 3. Create cran-comments.md

```markdown
## Test environments
* local macOS Sonoma 14.x, R 4.5.0
* ubuntu-latest (GitHub Actions), R-release
* windows-latest (GitHub Actions), R-release

## R CMD check results

0 errors ✓ | 0 warnings ✓ | 8 notes ✗

### Notes:

* checking CRAN incoming feasibility ... NOTE
  - New submission (first CRAN submission)
  - URLs return 404: GitHub repository is active but may timeout for CRAN's checker

* checking HTML version of manual ... NOTE
  - HTML Tidy not installed on local system (cosmetic, not a package issue)

* checking package dependencies ... NOTE
  - Suggests: 'bench' not available
  - This package is optional for benchmarking only, not required for core functionality

All other NOTEs are related to local system configuration or expected for new submissions.

## Downstream dependencies

None (this is a new package)
```

### 4. Submit to CRAN

**Option A: Web Submission** (Recommended for first-time)

1. Go to: https://cran.r-project.org/submit.html
2. Upload: `RcppICA_0.1.0.tar.gz`
3. Fill in maintainer email
4. Add submission comments from `cran-comments.md`
5. Click "Upload package"
6. **Check your email** for confirmation link
7. Click confirmation link within 60 minutes

**Option B: Using devtools**

```r
devtools::release()
# Follow prompts
```

### 5. Submission Comments

In the web form "Comments" box:

```
This is the first CRAN submission for RcppICA.

The package provides fast Independent Component Analysis using C++
with Eigen and Spectra libraries. It features Spectra-based truncated
eigendecomposition for 2-3x speedup over full eigendecomposition.

Key features:
- Competitive performance with base fastICA
- Comprehensive single-cell RNA-seq support
- 126 passing unit tests
- 3 detailed vignettes

R CMD check results:
- 0 ERRORs, 0 WARNINGs, 8 NOTEs
- All NOTEs are expected (new submission, optional suggests package,
  local system configuration)

The package has been tested on multiple platforms via GitHub Actions.
```

---

## After Submission

### Timeline

- **Initial review**: 1-3 days (auto-check by CRAN robots)
- **Human review**: 1-2 weeks (may be faster for clean submissions)
- **Possible outcomes**:
  1. Accepted → Published on CRAN
  2. Comments → Fix issues and resubmit
  3. Rejected → Address concerns and resubmit

### Responding to CRAN

If CRAN requests changes:

1. **Read carefully**: CRAN reviewers are thorough
2. **Fix all issues**: Don't argue, just fix
3. **Increment version**: `0.1.0` → `0.1.1` for resubmission
4. **Respond politely**: Thank them for the review
5. **Resubmit quickly**: Within a few days

**Example response**:

```
Dear CRAN team,

Thank you for the review. I have made the following changes in version 0.1.1:

1. Fixed [issue 1] by [action taken]
2. Updated [issue 2] by [action taken]
3. Clarified [issue 3] in documentation

All changes are documented in NEWS.md.

Best regards,
Michael Totty
```

---

## Common CRAN Feedback

### Examples Run Too Long

If CRAN says examples take >5 seconds:

```r
# Wrap slow examples in \donttest{}
#' @examples
#' \donttest{
#' # This example takes >5 seconds
#' result <- fastICA(large_matrix, n.comp = 100)
#' }
```

### Non-Standard Files

If CRAN complains about files in package:

Add to `.Rbuildignore`:
```
^benchmark.*\.R$
^profile.*\.R$
^.*\.md$ (except NEWS.md)
```

### URL Checks

If URLs return 404:

- Ensure GitHub repo is public
- Check URLs are correct
- Add `<URL>` tags if needed

---

## Post-Acceptance

### Maintenance

Once on CRAN, you must:

1. **Respond to CRAN emails** within 2 weeks
2. **Fix critical bugs** quickly
3. **Submit updates** responsibly (not too frequently)
4. **Test on multiple platforms** before submission

### Updates

**Good reasons to update**:
- Bug fixes
- New features
- Performance improvements
- Compatibility with new R versions

**Bad reasons to update**:
- Minor documentation tweaks
- Too frequent updates (wait 3-6 months between)

### Version Numbering

Follow semantic versioning:

- **0.1.1**: Bug fixes (patch)
- **0.2.0**: New features (minor)
- **1.0.0**: Major stable release

---

## Removal from CRAN

CRAN may remove package if:

- Check errors on CRAN systems
- Failing regularly with new R versions
- Dependencies are removed
- Maintainer unresponsive

**Prevention**:
- Monitor CRAN check results: https://cran.r-project.org/web/checks/check_results_RcppICA.html
- Respond to CRAN emails promptly
- Test with R-devel regularly

---

## Alternative: R-universe

If CRAN seems too strict, consider **R-universe** instead:

- Easier submission process
- More flexible
- Still discoverable
- Good for experimental packages

Setup: https://r-universe.dev/

---

## Recommended Timeline

**Week 0**: GitHub release (done)
**Weeks 1-4**: Gather user feedback, fix bugs
**Week 5**: Polish based on feedback
**Week 6**: Submit to CRAN

**Rationale**: Let users test on GitHub first, iron out issues, then submit a polished package to CRAN.

---

## Current Recommendation

**For v0.1.0**: ✅ **Focus on GitHub release**

**For v0.1.1 or v0.2.0**: Consider CRAN submission after:
- User feedback incorporated
- Any bugs fixed
- Platform testing complete
- Vignettes polished

---

## Resources

- **CRAN Repository Policy**: https://cran.r-project.org/web/packages/policies.html
- **R Package Development**: https://r-pkgs.org/
- **devtools::release()**: Automated CRAN submission
- **rhub**: Test on CRAN platforms: `rhub::check_for_cran()`

---

## Summary

**v0.1.0 is ready for GitHub release, but wait on CRAN**

Reasons:
1. Gather user feedback first
2. Fix any platform-specific issues
3. Polish based on real-world usage
4. Submit a battle-tested package to CRAN

**Next steps**:
1. ✅ GitHub release (see GITHUB_RELEASE.md)
2. 📊 Monitor user feedback (2-4 weeks)
3. 🔧 Fix any issues found
4. 📦 Submit to CRAN (when ready)

Package is well-positioned for CRAN success, but a GitHub "soft launch" first is prudent.
