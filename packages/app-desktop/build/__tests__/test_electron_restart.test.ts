import { afterEach, describe, expect, test } from "vitest";

import { restartChain, type ElectronChild } from "../electron-restart";

const pendingDeadlines: Array<() => void> = [];

const NEVER = () =>
  new Promise<void>((resolve) => {
    pendingDeadlines.push(resolve);
  });

const NOW = () => Promise.resolve();

describe("electron_restart", () => {
  afterEach(() => {
    pendingDeadlines.splice(0).forEach((resolve) => {
      resolve();
    });
  });

  test("a restart stops the running app, waits for its exit, then starts the next one", async () => {
    const log: string[] = [];
    const child = fakeChild(log, "SIGTERM");

    await restartChain()({ current: () => child, start: logged(log, "start"), exitDeadline: NEVER });

    expect(log).toStrictEqual(["kill SIGTERM", "exit", "start"]);
  });

  test("the exit listener of the plugin is removed, so the stop does not end the dev server", async () => {
    const log: string[] = [];
    const child = fakeChild(log, "SIGTERM");

    child.once("exit", () => {
      log.push("plugin exit listener");
    });
    await restartChain()({ current: () => child, start: logged(log, "start"), exitDeadline: NEVER });

    expect(log).toStrictEqual(["kill SIGTERM", "exit", "start"]);
  });

  test("an app that does not exit before the deadline is killed with SIGKILL", async () => {
    const log: string[] = [];
    const child = fakeChild(log, "SIGKILL");

    await restartChain()({ current: () => child, start: logged(log, "start"), exitDeadline: NOW });

    expect(log).toStrictEqual(["kill SIGTERM", "kill SIGKILL", "exit", "start"]);
  });

  test("a first start with no running app only starts", async () => {
    const log: string[] = [];

    await restartChain()({ current: () => undefined, start: logged(log, "start") });

    expect(log).toStrictEqual(["start"]);
  });

  test("two quick restarts run one after another", async () => {
    const log: string[] = [];
    const restart = restartChain();
    let child: ElectronChild | undefined;

    const start = (name: string) => async () => {
      log.push(`start ${name}`);
      await Promise.resolve();
      child = fakeChild(log, "SIGTERM");
    };

    await Promise.all([
      restart({ current: () => child, start: start("first"), exitDeadline: NEVER }),
      restart({ current: () => child, start: start("second"), exitDeadline: NEVER }),
    ]);

    expect(log).toStrictEqual(["start first", "kill SIGTERM", "exit", "start second"]);
  });
});

function fakeChild(log: string[], exitsOn: NodeJS.Signals): ElectronChild {
  let listeners: Array<() => void> = [];
  let exitCode: number | null = null;

  return {
    get exitCode() {
      return exitCode;
    },
    signalCode: null,
    removeAllListeners() {
      listeners = [];
    },
    once(_event, listener) {
      listeners.push(listener);
    },
    kill(signal = "SIGTERM") {
      log.push(`kill ${signal}`);

      if (signal === exitsOn) {
        exitCode = 0;
        log.push("exit");
        listeners.forEach((listener) => {
          listener();
        });
      }

      return true;
    },
  };
}

function logged(log: string[], entry: string): () => Promise<void> {
  return async () => {
    log.push(entry);
    await Promise.resolve();
  };
}
