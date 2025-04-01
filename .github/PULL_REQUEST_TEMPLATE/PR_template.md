### What's new
- Describe what was changed in the code and why it is useful or necessary
- Closes #<issue_number>

```[tasklist]
### Blocked by...
- [ ] merge after user/repo_name#PR_number
- [ ] depends on #PR_number
```

### Run it locally 
1. Checkout the PR branch.

> Assuming 1) you already have a local copy of `q2-annotate` that is installed in editable mode in your activated virtual environment and 2) the working directory is `q2-annotate`; run the following. 

```bash
PR=<PR_number>
git fetch origin pull/${PR}/head:pr-${PR}
git checkout pr-${PR}
```

2. Let's get you some data to play with: 
```bash
<your code here>
```

3. Test it out!
```bash
<your code here>
```

### TODO
- [ ] Fill in *Whats new* section. Erase any bullet points that are not used.
- [ ] Erase *Blocked by...* task list if not used.
- [ ] In *Run it locally*, **step 1**, make sure to write the PR number after you have submitted the PR.
- [ ] In *Run it locally*, **step 2**, fill in the code to fetch the inputs to test the code.
- [ ] In *Run it locally*, **step 3**, fill in the code to test the changes. 
- [ ] Erase the *TODO* section.
