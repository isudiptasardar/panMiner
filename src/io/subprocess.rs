//! Centralized subprocess execution with timeout support.
//!
//! All external tool invocations should use `run_with_timeout()` to prevent
//! hung subprocesses from blocking the pipeline indefinitely.

use std::io::Read;
use std::process::{Command, Output, Stdio};
use std::time::{Duration, Instant};

use crate::error::{Error, Result};

/// Default subprocess timeout in seconds (1 hour).
const DEFAULT_TIMEOUT_SECS: u64 = 3600;

/// Run a subprocess command with a configurable timeout.
///
/// If the subprocess does not exit within `timeout`, it is killed and
/// `Error::SubprocessTimeout` is returned.
///
/// # Arguments
///
/// * `command` - The `Command` to execute (already configured with args, env, etc.)
/// * `tool_name` - Human-readable tool name for error messages
/// * `timeout_secs` - Optional custom timeout in seconds (default: 3600)
///
/// # Returns
///
/// The `Output` on success, or an error if the process fails or times out.
pub fn run_with_timeout(
    command: &mut Command,
    tool_name: &str,
    timeout_secs: Option<u64>,
) -> Result<Output> {
    let timeout = Duration::from_secs(timeout_secs.unwrap_or(DEFAULT_TIMEOUT_SECS));
    let tool = tool_name.to_string();

    command.stdout(Stdio::piped());
    command.stderr(Stdio::piped());

    let start = Instant::now();
    let mut child = command.spawn().map_err(|e| {
        Error::ExternalTool(format!("Failed to spawn {}: {}", tool, e))
    })?;

    loop {
        match child.try_wait() {
            Ok(Some(status)) => {
                // Process exited — read remaining stdout/stderr
                let stdout = match child.stdout.take() {
                    Some(mut h) => {
                        let mut buf = Vec::new();
                        let _ = h.read_to_end(&mut buf);
                        buf
                    }
                    None => Vec::new(),
                };
                let stderr = match child.stderr.take() {
                    Some(mut h) => {
                        let mut buf = Vec::new();
                        let _ = h.read_to_end(&mut buf);
                        buf
                    }
                    None => Vec::new(),
                };
                return Ok(Output { status, stdout, stderr });
            }
            Ok(None) => {
                if start.elapsed() >= timeout {
                    let _ = child.kill();
                    let _ = child.wait(); // reap the zombie
                    return Err(Error::SubprocessTimeout {
                        tool: tool.clone(),
                        timeout_secs: timeout.as_secs(),
                    });
                }
                std::thread::sleep(Duration::from_millis(100));
            }
            Err(e) => {
                return Err(Error::ExternalTool(format!(
                    "Error waiting for {}: {}",
                    tool, e
                )));
            }
        }
    }
}

/// Get the default timeout duration.
#[allow(dead_code)]
pub fn default_timeout() -> Duration {
    Duration::from_secs(DEFAULT_TIMEOUT_SECS)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_timeout_value() {
        assert_eq!(DEFAULT_TIMEOUT_SECS, 3600);
        assert_eq!(default_timeout().as_secs(), 3600);
    }

    #[test]
    fn test_run_with_timeout_success() {
        let mut cmd = Command::new("echo");
        cmd.arg("hello");
        let result = run_with_timeout(&mut cmd, "echo", Some(10));
        assert!(result.is_ok());
        let output = result.unwrap();
        assert!(output.status.success());
    }

    #[test]
    fn test_run_with_timeout_nonexistent_command() {
        let mut cmd = Command::new("nonexistent_tool_xyz_12345");
        let result = run_with_timeout(&mut cmd, "nonexistent_tool_xyz_12345", Some(5));
        assert!(result.is_err());
    }

    #[test]
    fn test_run_with_timeout_custom_timeout() {
        let mut cmd = Command::new("echo");
        cmd.arg("fast");
        let result = run_with_timeout(&mut cmd, "echo", Some(1));
        assert!(result.is_ok());
    }

    #[test]
    fn test_run_with_timeout_expires() {
        // Use a long-running command with a 1-second timeout
        #[cfg(target_os = "windows")]
        let mut cmd = Command::new("ping");
        #[cfg(target_os = "windows")]
        cmd.args(["-n", "60", "127.0.0.1"]);

        #[cfg(not(target_os = "windows"))]
        let mut cmd = Command::new("sleep");
        #[cfg(not(target_os = "windows"))]
        cmd.arg("60");

        let result = run_with_timeout(&mut cmd, "slow_cmd", Some(1));
        assert!(result.is_err());
        match result.unwrap_err() {
            Error::SubprocessTimeout { tool, timeout_secs } => {
                assert_eq!(tool, "slow_cmd");
                assert_eq!(timeout_secs, 1);
            }
            other => panic!("Expected SubprocessTimeout, got: {}", other),
        }
    }
}