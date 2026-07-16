#!/usr/bin/env python3
import os
import re
import sys

# Directory containing markdown files
dir_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

invalid_left_right_pattern = re.compile(r'\\(?:left|right)\s*(?<!\\)[{}]')

def check_braces_balance(equation):
    stripped = re.sub(r'\\{', '', equation)
    stripped = re.sub(r'\\}', '', stripped)
    open_count = stripped.count('{')
    close_count = stripped.count('}')
    return open_count == close_count, open_count, close_count

def check_environments(equation):
    begins = re.findall(r'\\begin\{([a-zA-Z*]+)\}', equation)
    ends = re.findall(r'\\end\{([a-zA-Z*]+)\}', equation)
    if len(begins) != len(ends):
        return False, f"Mismatched environment count: begins={begins}, ends={ends}"
    
    stack = []
    for token in re.findall(r'\\(?:begin|end)\{[a-zA-Z*]+\}', equation):
        if 'begin' in token:
            env = re.search(r'\\begin\{([a-zA-Z*]+)\}', token).group(1)
            stack.append(env)
        else:
            env = re.search(r'\\end\{([a-zA-Z*]+)\}', token).group(1)
            if not stack:
                return False, f"Unexpected \\end{{{env}}} without matching \\begin"
            last = stack.pop()
            if last != env:
                return False, f"Mismatched environments: \\begin{{{last}}} closed by \\end{{{env}}}"
    return True, ""

def validate_equation(equation):
    errors = []
    invalid_lr = invalid_left_right_pattern.findall(equation)
    if invalid_lr:
        errors.append(f"Invalid left/right delimiter: found raw braces {invalid_lr} (must use \\{{ and \\}})")
    balanced, opens, closes = check_braces_balance(equation)
    if not balanced:
        errors.append(f"Unbalanced braces: {opens} open braces vs {closes} close braces")
    env_ok, env_msg = check_environments(equation)
    if not env_ok:
        errors.append(env_msg)
    if "\\operatorname" in equation:
        errors.append("Unsupported \\operatorname found (use \\mathrm or \\text instead)")
    return errors

def scan_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
        
    errors_found = []
    
    # 1. Block math
    block_pattern = re.compile(r'```math\n(.*?)\n```', re.DOTALL)
    for match in block_pattern.finditer(content):
        equation = match.group(1)
        start_char = match.start()
        line_num = content[:start_char].count('\n') + 1
        errs = validate_equation(equation)
        if errs:
            errors_found.append((line_num, "Block", equation, errs))
            
    # 2. Inline math
    inline_pattern = re.compile(r'(?<!\\)\$(?:`\s*(.*?)\s*`|(?!`)(.*?))(?<!\\)\$', re.DOTALL)
    for match in inline_pattern.finditer(content):
        start_char = match.start()
        before = content[:start_char]
        if before.count('```') % 2 != 0:
            continue
            
        equation = match.group(1) or match.group(2)
        if not equation or equation.strip() == "":
            continue
            
        line_num = before.count('\n') + 1
        errs = validate_equation(equation)
        if errs:
            errors_found.append((line_num, "Inline", equation, errs))
            
    return errors_found

def main():
    total_files = 0
    total_errors = 0
    
    for root, dirs, files in os.walk(dir_path):
        # Exclude hidden directories like .git
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        
        for file in files:
            if file.endswith('.md'):
                total_files += 1
                filepath = os.path.join(root, file)
                errs = scan_file(filepath)
                if errs:
                    total_errors += len(errs)
                    rel_path = os.path.relpath(filepath, dir_path)
                    print(f"\n[ERROR] in file: {rel_path}")
                    for line_num, eq_type, eq, problems in errs:
                        print(f"  Line {line_num} ({eq_type} math):")
                        for p in problems:
                            print(f"    - {p}")
                        print(f"    Code: {eq.strip()}")
                        
    if total_errors > 0:
        print(f"\nValidation failed with {total_errors} error(s) in {total_files} files.")
        sys.exit(1)
    else:
        print(f"\nValidation successful: Checked {total_files} markdown files. No errors found!")
        sys.exit(0)

if __name__ == '__main__':
    main()
